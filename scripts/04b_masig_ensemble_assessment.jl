using BlackBoxOptim
using PairPlots

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "_masig_config.jl"))

reef_id = "masig"
ensemble_dir = joinpath(OUTPUT_DIR, "ensemble", "torres_strait", "masig")

ensemble_output = deserialize(joinpath(ensemble_dir, "$(reef_id)_ensemble_output.dat"))
masig_ensemble = deserialize(joinpath(ensemble_dir, "$(reef_id)_tracked_candidates.dat"))
ensemble_params = hcat(masig_ensemble.candidates...)

parameter_identifiability_metrics(ensemble_params, ENSEMBLE_PARAM_NAMES)

corr_df = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.5
)
corr_df[!, :Correlation] .= round.(corr_df.Correlation; digits=3)
CSV.write(joinpath(ensemble_dir, "$(reef_id)_parameter_correlations.csv"), corr_df)

# ── Create paths ──────────────────────────────────────────────────────────────
ensemble_data_dir = joinpath(OUTPUT_DIR, "sensitivity", "torres_strait", "masig", "ensemble")
mkpath(ensemble_data_dir)

fn_unconstrained_samples = joinpath(ensemble_data_dir, "$(reef_id)_unconstrained_samples.dat")
fn_unconstrained_fitness  = joinpath(ensemble_data_dir, "$(reef_id)_unconstrained_fitness.dat")
fn_unconstrained_pawn     = joinpath(ensemble_data_dir, "$(reef_id)_unconstrained_pawn_results.dat")

ensemble_fig_dir = joinpath(FIG_DIR, "sensitivity", "torres_strait", "masig", "ensemble")
mkpath(ensemble_fig_dir)

# Load observations and rebuild sim_year_range
reef_df = CSV.read(joinpath(OUTPUT_DIR, "Masig_Reef_EcoRRAP_estimate.csv"), DataFrame)
reef_obs = DataFrame(;
    SAMPLE_DATE=reef_df.report_year,
    MEAN_LIVE_CORAL=reef_df.mean,
    LOWER=reef_df.lower,
    UPPER=reef_df.upper
)

benthic_estimate = CSV.read(joinpath(OUTPUT_DIR, "ecorrap_benthic", "masig_estimate.csv"), DataFrame)
sim_start_year = benthic_estimate.year[1]
sim_end_year = benthic_estimate.year[end]

sim_year_range = create_simulation_dates(
    sim_start_year, sim_end_year, calib_settings.start_month
)

reef_state = ensemble_output.reef_state
env_conditions = ensemble_output.env_conditions

sim_indices, ref_indices, matched_dates = find_closest_dates(
    collect(Date.(sim_year_range)),
    Date.(reef_obs.SAMPLE_DATE);
    max_days=calib_settings.date_match_max_days
)

sim_indices, ref_indices, matched_dates = exclude_years_from_indices(
    sim_indices, ref_indices, matched_dates, reef_config.exclude_years
)

year_span = year.(matched_dates)
sim_benthic_years = [year_span .∈ Ref(benthic_estimate.year)][1]
aligned_years = year_span[sim_benthic_years]
benthic_data = benthic_estimate[benthic_estimate.year .∈ Ref(aligned_years), :]

sample_method = SobolSample(; R=OwenScramble(; base=2, pad=32))
n = 8192  # 2^13

n_threads = Threads.nthreads()
runner_pool = create_objective_pool(
    reef_state, env_conditions, reef_obs.MEAN_LIVE_CORAL,
    sim_indices, ref_indices, reef_config.area,
    sim_benthic_years, benthic_data,
    opt_config.random_seed, calib_settings.use_scalers
)

# ── Unconstrained sensitivity analysis ───────────────────────────────────────
if !isfile(fn_unconstrained_samples)
    unif_dists = [Uniform(l, u) for (l, u) in param_bounds]
    n_bounds = length(param_bounds)
    unc_samples = Matrix(
        QMC.sample(n, zeros(n_bounds), ones(n_bounds), sample_method)'
    )

    unc_samples[:, 2:6] = mapslices(gamma_to_dirichlet, unc_samples[:, 2:6]; dims=2)

    free_vary_cols = [1, 7:23...]
    unc_samples[:, free_vary_cols] = Matrix(
        Distributions.quantile.(unif_dists[[1, 7:23...]], unc_samples[:, free_vary_cols]')'
    )

    serialize(fn_unconstrained_samples, unc_samples)

    @info "Running unconstrained sample ($(n_threads) threads)"
    unc_fitness_scores = Vector{Float64}(undef, size(unc_samples, 1))
    Threads.@threads :static for i in axes(unc_samples, 1)
        unc_fitness_scores[i] = runner_pool[Threads.threadid()](collect(unc_samples[i, :]))
    end

    serialize(fn_unconstrained_fitness, unc_fitness_scores)

    unc_pawn_sa_results = pawn(unc_samples, unc_fitness_scores, ENSEMBLE_PARAM_NAMES)
    serialize(fn_unconstrained_pawn, unc_pawn_sa_results)

    f = plot_pawn_heatmap(
        unc_pawn_sa_results, "Unconstrained SA - $(reef_id)";
        stats=[:mean, :std], xticklabelrotation=π / 2
    )
    save(joinpath(ensemble_fig_dir, "$(reef_id)_unconstrained_sa.png"), f; px_per_unit=DPI)
else
    unc_samples = deserialize(fn_unconstrained_samples)
    unc_fitness_scores = deserialize(fn_unconstrained_fitness)
    unc_pawn_sa_results = deserialize(fn_unconstrained_pawn)
end

# ── Constrained sensitivity analysis ─────────────────────────────────────────
fn_constrained_samples = joinpath(ensemble_data_dir, "$(reef_id)_constrained_samples.dat")
fn_constrained_fitness  = joinpath(ensemble_data_dir, "$(reef_id)_constrained_fitness.dat")
fn_constrained_pawn     = joinpath(ensemble_data_dir, "$(reef_id)_constrained_pawn_results.dat")

if !isfile(fn_constrained_samples)
    cons_param_bounds = extrema.(eachrow(ensemble_params))
    unif_dists = [Uniform(l, u) for (l, u) in cons_param_bounds]
    n_bounds = length(cons_param_bounds)
    cons_samples = Matrix(
        QMC.sample(n, zeros(n_bounds), ones(n_bounds), sample_method)'
    )

    cons_samples[:, 2:6] = mapslices(gamma_to_dirichlet, cons_samples[:, 2:6]; dims=2)

    free_vary_cols = [1, 7:23...]
    cons_samples[:, free_vary_cols] = Matrix(
        Distributions.quantile.(unif_dists[[1, 7:23...]], cons_samples[:, free_vary_cols]')'
    )
    serialize(fn_constrained_samples, cons_samples)

    @info "Running ensemble-constrained sample ($(n_threads) threads)"
    cons_fitness_scores = Vector{Float64}(undef, size(cons_samples, 1))
    Threads.@threads :static for i in axes(cons_samples, 1)
        cons_fitness_scores[i] = runner_pool[Threads.threadid()](collect(cons_samples[i, :]))
    end

    serialize(fn_constrained_fitness, cons_fitness_scores)

    cons_pawn_sa_results = pawn(cons_samples, cons_fitness_scores, ENSEMBLE_PARAM_NAMES)
    serialize(fn_constrained_pawn, cons_pawn_sa_results)

    f = plot_pawn_heatmap(
        cons_pawn_sa_results, "Constrained SA - $(reef_id)";
        stats=[:mean, :std], xticklabelrotation=π / 2
    )
    save(joinpath(ensemble_fig_dir, "$(reef_id)_constrained_sa.png"), f; px_per_unit=DPI)
else
    cons_samples = deserialize(fn_constrained_samples)
    cons_fitness_scores = deserialize(fn_constrained_fitness)
    cons_pawn_sa_results = deserialize(fn_constrained_pawn)
end

# ── Export PAWN results as CSV ────────────────────────────────────────────────
function _pawn_to_csv(results::YAXArray, path::String)
    factors = string.(collect(results.axes[1]))
    stats   = collect(results.axes[2])
    df = DataFrame(:parameter => factors)
    for s in stats
        df[!, string(s)] = collect(results[PAWNᵢ=At(s)])
    end
    CSV.write(path, df)
    return df
end

_pawn_to_csv(
    unc_pawn_sa_results,
    joinpath(ensemble_data_dir, "$(reef_id)_unconstrained_pawn_results.csv")
)
_pawn_to_csv(
    cons_pawn_sa_results,
    joinpath(ensemble_data_dir, "$(reef_id)_constrained_pawn_results.csv")
)

# ── Publication theme ─────────────────────────────────────────────────────────
fontsize_theme = Theme(; fontsize=14)
set_theme!(fontsize_theme)

# ── Ensemble correlations ─────────────────────────────────────────────────────
parameter_identifiability_metrics(ensemble_params, ENSEMBLE_PARAM_NAMES)
corr_threshold = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.5
)

corr_params = unique(vcat(corr_threshold.Param1, corr_threshold.Param2))

df = DataFrame(Matrix(ensemble_params'), :auto)
rename!(df, ENSEMBLE_PARAM_NAMES)

target_df = df[:, corr_params]
f = pairplot(
    target_df => (
        PairPlots.Series(target_df),
        PairPlots.HexBin(; colormap=Makie.cgrad([:transparent, :black])),
        PairPlots.Scatter(; alpha=0.5),
        PairPlots.Contour(),
        PairPlots.MarginDensity(; bandwidth=0.1),
        PairPlots.MarginQuantileText()
    )
)

resize_to_layout!(f)
sleep(5)
save(joinpath(ensemble_fig_dir, "$(reef_id)_ensemble_corr_param_pairplot.png"), f; px_per_unit=DPI)

# ── Identify most influential parameters from constrained PAWN ────────────────
most_influential = collect(
    sortperm(cons_pawn_sa_results[PAWNᵢ=At(:mean)]; rev=true)[1:10]
)

factor_names = collect(collect(cons_pawn_sa_results.factors[most_influential]))

target_df = df[:, factor_names]

# ── Uncertainty reduction: prior vs. posterior ────────────────────────────────
posterior_samples = Matrix(ensemble_params')

σ_prior     = vec(std(unc_samples;       dims=1))
σ_posterior = vec(std(posterior_samples; dims=1))
rel_reduction = (σ_prior .- σ_posterior) ./ σ_prior

# Top-3 most-sensitive (PAWN) parameters
top3_sensitive = most_influential[1:3]

function plot_uncertainty_panel!(ax, pidx, unc_samples, posterior_samples, rel_reduction)
    unc_v  = unc_samples[:,       pidx]
    post_v = posterior_samples[:, pidx]
    pct    = round(rel_reduction[pidx] * 100; digits=1)

    hist!(ax, unc_v;  normalization=:pdf, bins=50, color=(:steelblue, 0.45))
    hist!(ax, post_v; normalization=:pdf, bins=50, color=(:orangered, 0.45))
    density!(ax, unc_v;  color=(:steelblue, 0.0), strokecolor=:steelblue, strokewidth=2.5)
    density!(ax, post_v; color=(:orangered, 0.0), strokecolor=:orangered, strokewidth=2.5)

    σ_prior_v = round(std(unc_v);  digits=3)
    σ_post_v  = round(std(post_v); digits=3)
    text!(ax, 0.97, 0.97;
        text="σ prior: $(σ_prior_v) → $(σ_post_v) (−$(pct)%)",
        align=(:right, :top),
        space=:relative,
        fontsize=18,
        color=:black
    )
end

# ── Combined figure: constrained SA pairplot (a) + uncertainty reduction (b) ──
fig_combined = Figure(; size=(2000, 2750))

pairplot(
    fig_combined[1, 1],
    target_df => (
        PairPlots.Series(target_df),
        PairPlots.HexBin(; colormap=Makie.cgrad([:transparent, :black])),
        PairPlots.Scatter(; alpha=0.5),
        PairPlots.Contour(),
        PairPlots.MarginDensity(; bandwidth=0.1),
        PairPlots.MarginQuantileText()
    );
    labels=Dict(fn => rich(string(fn); fontsize=20) for fn in factor_names)
)

gl_unc = fig_combined[2, 1] = GridLayout()
Label(gl_unc[1, 0]; text="Most\ninfluential\n(PAWN)", tellheight=false, rotation=π/2, fontsize=22)

for (col, pidx) in enumerate(top3_sensitive)
    local ax = Axis(
        gl_unc[1, col];
        xlabel=string(ENSEMBLE_PARAM_NAMES[pidx]),
        ylabel=col == 1 ? "Density" : "",
        xlabelsize=22,
        ylabelsize=22,
        xticklabelsize=20,
        yticklabelsize=20,
    )
    plot_uncertainty_panel!(ax, pidx, unc_samples, posterior_samples, rel_reduction)
end

Legend(gl_unc[0, 1:3],
    [PolyElement(; color=(:steelblue, 0.6)), PolyElement(; color=(:orangered, 0.6))],
    ["Prior (unconstrained)", "Posterior (calibrated ensemble)"];
    orientation=:horizontal,
    framevisible=false
)
Label(gl_unc[0, 0]; text="", tellheight=false)

Label(fig_combined[1, 1, TopLeft()], "(A)"; fontsize=16, font=:bold, padding=(4, 0, 4, 0))
Label(fig_combined[2, 1, TopLeft()], "(B)"; fontsize=16, font=:bold, padding=(4, 0, 4, 0))

rowsize!(fig_combined.layout, 1, 1950)
rowsize!(fig_combined.layout, 2, 350)

sleep(5)
save(
    joinpath(ensemble_fig_dir, "$(reef_id)_cons_ensemble_sa_param_pairplot.png"),
    fig_combined; px_per_unit=DPI
)

# ═══════════════════════════════════════════════════════════════════════════════
# ADDITION 2 — Temporal (time-varying) PAWN sensitivity analysis
# ═══════════════════════════════════════════════════════════════════════════════

"""
    create_timeseries_function(reef_state, env_conditions, area, seed, use_scalers)

Create a function that runs the model and returns the full simulated coral cover
trajectory rather than aggregating to a scalar fitness score.  Thread-safe via
an internal pool of reef_state copies.
"""
function create_timeseries_function(
    reef_state::ReefState,
    env_conditions::YAXArray,
    area::Float32,
    seed::Int,
    use_scalers::Bool
)
    n_threads = Threads.nthreads()
    state_pool = [copy(reef_state) for _ in 1:n_threads]
    rng_pool = [
        Random.seed!(copy(Random.default_rng()), seed + i) for i in 0:(n_threads - 1)
    ]

    function timeseries(x; return_by_group::Bool=false)
        tid = Threads.threadid()
        rs = state_pool[tid]
        rng = rng_pool[tid]

        x = vcat(x[1], gamma_to_dirichlet(x[2:6]), x[7:end])

        CoralFlow.set_population!(rs, x)

        if use_scalers
            n_grps = CoralFlow.n_groups(rs)
            scaler_end = 17 + n_grps - 1
            loc_scalers = x[17:scaler_end]
            CoralFlow.assign_scalers!(rs, loc_scalers)
            recruitment_proportion = x[end - 1]
            self_seeding_proportion = x[end]
            CoralFlow.run_model!(
                rs, env_conditions;
                recruits=Float32(recruitment_proportion),
                self_seed=Float32(self_seeding_proportion),
                rng=rng
            )
        else
            CoralFlow.run_model!(rs, env_conditions; rng=rng)
        end

        if return_by_group
            return CoralFlow.group_cover_timeseries(rs) ./ area
        else
            return CoralFlow.coral_cover(rs) ./ area
        end
    end

    return timeseries
end

# ── Output paths ──────────────────────────────────────────────────────────────
fn_temporal_ts   = joinpath(ensemble_data_dir, "$(reef_id)_temporal_ts_outputs.dat")
fn_temporal_pawn = joinpath(ensemble_data_dir, "$(reef_id)_temporal_pawn_results.dat")

# ── Recover per-timestep trajectories ─────────────────────────────────────────
if !isfile(fn_temporal_ts)
    if hasproperty(ensemble_output, :trajectories) && !isnothing(ensemble_output.trajectories)
        @info "Loading trajectories from ensemble_output"
        ts_outputs =
            ensemble_output.trajectories isa Matrix ?
            ensemble_output.trajectories :
            Matrix(ensemble_output.trajectories')
    else
        @info "Re-running model on $(size(ensemble_params, 2)) ensemble candidates"

        n_candidates = size(ensemble_params, 2)
        n_sim_steps = length(env_conditions.timestep)

        if n_candidates < 100
            @warn "Only $(n_candidates) ensemble candidates available. " *
                "Temporal PAWN estimates may be unreliable."
        end

        model_ts_runner = create_timeseries_function(
            reef_state, env_conditions, reef_config.area,
            opt_config.random_seed, calib_settings.use_scalers
        )

        ts_outputs = Matrix{Float32}(undef, n_candidates, n_sim_steps)
        Threads.@threads :static for i in 1:n_candidates
            ts_outputs[i, :] = model_ts_runner(ensemble_params[:, i])
        end
    end

    serialize(fn_temporal_ts, ts_outputs)
else
    ts_outputs = deserialize(fn_temporal_ts)
end

n_sim_steps = size(ts_outputs, 2)
n_params = size(ensemble_params, 1)

X = Matrix(ensemble_params')

# ── Compute PAWN at each timestep ─────────────────────────────────────────────
if !isfile(fn_temporal_pawn)
    temporal_pawn = Matrix{Float64}(undef, n_params, n_sim_steps)

    @info "Computing temporal PAWN ($(n_params) params × $(n_sim_steps) timesteps)"
    for t in 1:n_sim_steps
        pawn_t = pawn(X, ts_outputs[:, t], ENSEMBLE_PARAM_NAMES)
        temporal_pawn[:, t] = pawn_t[PAWNᵢ=At(:mean)].data
    end

    serialize(fn_temporal_pawn, temporal_pawn)
else
    temporal_pawn = deserialize(fn_temporal_pawn)
end

decimal_years = collect(sim_start_year:sim_end_year)
disturbance_years = Float64.(reef_config.disturbance_years)

# ── Figure 1: heatmap (parameters × time) ────────────────────────────────────
row_order = sortperm(vec(mean(temporal_pawn; dims=2)); rev=true)
sorted_pawn = temporal_pawn[row_order, :]
sorted_names = ENSEMBLE_PARAM_NAMES[row_order]

fig_tsa = Figure(; size=(1100, max(300, 28 * n_params)))
ax_tsa = Axis(
    fig_tsa[1, 1];
    xlabel="Year",
    ylabel="Parameter",
    title="Temporal PAWN sensitivity — $(reef_id)",
    yticks=(1:n_params, string.(sorted_names)),
    yreversed=false
)

hm = heatmap!(
    ax_tsa, decimal_years, 1:n_params, sorted_pawn';
    colormap=:viridis,
    colorrange=(0.0, max(0.1, maximum(sorted_pawn)))
)

for yr in disturbance_years
    vlines!(ax_tsa, yr; color=(:red, 0.7), linewidth=1.5, linestyle=:dash)
end

Colorbar(fig_tsa[1, 2], hm; label="PAWN mean index", width=14)

save(
    joinpath(ensemble_fig_dir, "$(reef_id)_temporal_pawn_heatmap.png"),
    fig_tsa; px_per_unit=DPI
)

# ── Figure 2: line plot for the top-k parameters ──────────────────────────────
top_k = min(10, n_params)
top_k_rows = row_order[1:top_k]
top_k_names = ENSEMBLE_PARAM_NAMES[top_k_rows]

fig_lines = Figure(; size=(900, 400))
ax_lines = Axis(
    fig_lines[1, 1];
    xlabel="Year",
    ylabel="PAWN mean index",
    title="Temporal sensitivity — top $(top_k) parameters"
)

palette = Makie.wong_colors()
for (k, (pidx, pname)) in enumerate(zip(top_k_rows, top_k_names))
    lines!(ax_lines, decimal_years, temporal_pawn[pidx, :];
        label=string(pname),
        color=palette[mod1(k, length(palette))],
        linewidth=1.8
    )
end

y_max = maximum(temporal_pawn[top_k_rows, :])
for yr in disturbance_years
    vlines!(ax_lines, yr; color=(:red, 0.55), linewidth=1.2, linestyle=:dash)
    text!(ax_lines, yr + 0.1, y_max * 0.95;
        text="disturbance",
        rotation=π / 2,
        fontsize=9,
        color=(:red, 0.7)
    )
end

Legend(fig_lines[1, 2], ax_lines; framevisible=false)

save(
    joinpath(ensemble_fig_dir, "$(reef_id)_temporal_pawn_top$(top_k)_lines.png"),
    fig_lines; px_per_unit=DPI
)

# ── Per-group temporal PAWN ───────────────────────────────────────────────────
fn_temporal_ts_bygroup = joinpath(
    ensemble_data_dir, "$(reef_id)_temporal_ts_outputs_bygroup.dat"
)

if !isfile(fn_temporal_ts_bygroup)
    n_candidates = size(ensemble_params, 2)
    n_grps = CoralFlow.n_groups(reef_state)

    model_ts_runner = create_timeseries_function(
        reef_state, env_conditions, reef_config.area, opt_config.random_seed, true
    )

    ts_outputs_bygroup = Array{Float32,3}(undef, n_candidates, n_sim_steps, n_grps)
    Threads.@threads :static for i in 1:n_candidates
        ts_outputs_bygroup[i, :, :] = model_ts_runner(
            ensemble_params[:, i]; return_by_group=true
        )
    end
    serialize(fn_temporal_ts_bygroup, ts_outputs_bygroup)
else
    ts_outputs_bygroup = deserialize(fn_temporal_ts_bygroup)
end

# ═══════════════════════════════════════════════════════════════════════════════
# ADDITION 3 — Lagged PAWN sensitivity analysis
#
# For each lag k and window-start t, PAWN is computed on Y[i, t+k] directly:
# the level of cover k years after the window start.
# ═══════════════════════════════════════════════════════════════════════════════

lags_yr_all = [1, 2, 3, 5]
steps_per_yr = n_sim_steps / length(decimal_years)
lags_ts_all = round.(Int, lags_yr_all .* steps_per_yr)

# Drop lags that meet or exceed the simulation length — no valid windows would exist
valid_lag_mask = lags_ts_all .< n_sim_steps
lags_yr = lags_yr_all[valid_lag_mask]
lags_ts = lags_ts_all[valid_lag_mask]

if length(lags_yr) < length(lags_yr_all)
    @warn "Simulation length ($(n_sim_steps) steps) too short for lags $(lags_yr_all[.!valid_lag_mask]) yr — skipped"
end

fn_lagged_pawn = joinpath(ensemble_data_dir, "$(reef_id)_lagged_pawn_results.dat")

if !isfile(fn_lagged_pawn)
    lagged_pawn_results = Dict{Int,Matrix{Float64}}()

    for (k, lag_yr) in zip(lags_ts, lags_yr)
        n_valid = n_sim_steps - k
        lag_pawn = Matrix{Float64}(undef, n_params, n_valid)

        @info "Computing lagged PAWN — $(lag_yr)-yr window ($(k) timesteps)"
        for t in 1:n_valid
            pawn_t = pawn(X, ts_outputs[:, t + k], ENSEMBLE_PARAM_NAMES)
            lag_pawn[:, t] = pawn_t[PAWNᵢ=At(:mean)].data
        end

        lagged_pawn_results[k] = lag_pawn
    end

    serialize(fn_lagged_pawn, lagged_pawn_results)
else
    lagged_pawn_results = deserialize(fn_lagged_pawn)
end

# ── Figure per lag: heatmap (parameters × window-start year) ─────────────────
for (k, lag_yr) in zip(lags_ts, lags_yr)
    lag_pawn = lagged_pawn_results[k]
    n_valid = size(lag_pawn, 2)
    valid_years = decimal_years[1:n_valid]

    row_order_l = sortperm(vec(mean(lag_pawn; dims=2)); rev=true)
    sorted_pawn_l = lag_pawn[row_order_l, :]
    sorted_names_l = ENSEMBLE_PARAM_NAMES[row_order_l]

    fig_lag = Figure(; size=(1100, max(300, 28 * n_params)))
    ax_lag = Axis(
        fig_lag[1, 1];
        xlabel="Window start year",
        ylabel="Parameter",
        title="Temporal Sensitivity Analysis — cover at t+$(lag_yr)yr — $(reef_id)",
        yticks=(1:n_params, string.(sorted_names_l)),
        yreversed=false
    )

    hm_lag = heatmap!(
        ax_lag, valid_years, 1:n_params, sorted_pawn_l';
        colormap=:viridis,
        colorrange=(0.0, max(0.1, maximum(sorted_pawn_l)))
    )

    for yr in disturbance_years
        vlines!(ax_lag, yr; color=(:red, 0.7), linewidth=1.5, linestyle=:dash)
    end

    Colorbar(fig_lag[1, 2], hm_lag; label="PAWN mean index", width=14)

    save(
        joinpath(ensemble_fig_dir, "$(reef_id)_lagged_pawn_$(lag_yr)yr_heatmap.png"),
        fig_lag; px_per_unit=DPI
    )
end

# ── Cross-lag summary ──────────────────────────────────────────────────────────
fig_crosslag = Figure(; size=(900, 420))
ax_crosslag = Axis(
    fig_crosslag[1, 1];
    xlabel="Lag (years)",
    ylabel="Mean PAWN index",
    title="Parameter influence vs. recovery horizon — top $(top_k) — $(reef_id)",
    xticks=vcat(0, lags_yr)
)

for (rank, pidx) in enumerate(top_k_rows)
    lag0_mean = mean(temporal_pawn[pidx, :])
    lagged_means = [mean(lagged_pawn_results[k][pidx, :]) for k in lags_ts]
    all_lags = Float64.(vcat(0, lags_yr))
    all_means = vcat(lag0_mean, lagged_means)

    col = palette[mod1(rank, length(palette))]
    lines!(ax_crosslag, all_lags, all_means;
        label=string(ENSEMBLE_PARAM_NAMES[pidx]),
        color=col,
        linewidth=1.8
    )
    scatter!(ax_crosslag, all_lags, all_means; color=col, markersize=7)
end

Legend(fig_crosslag[1, 2], ax_crosslag; framevisible=false)

save(
    joinpath(ensemble_fig_dir, "$(reef_id)_lagged_pawn_crosslag_summary.png"),
    fig_crosslag; px_per_unit=DPI
)

# ── Per-group temporal PAWN heatmaps ──────────────────────────────────────────
group_names = CoralFlow.TARGET_GROUPS
for g in axes(ts_outputs_bygroup, 3)
    pawn_g = Matrix{Float64}(undef, n_params, n_sim_steps)
    for t in 1:n_sim_steps
        pawn_t = pawn(X, ts_outputs_bygroup[:, t, g], ENSEMBLE_PARAM_NAMES)
        pawn_g[:, t] = pawn_t[PAWNᵢ=At(:mean)].data
    end

    row_order_g = sortperm(vec(mean(pawn_g; dims=2)); rev=true)
    sorted_pawn_g = pawn_g[row_order_g, :]
    sorted_names_g = ENSEMBLE_PARAM_NAMES[row_order_g]

    fig_g = Figure(; size=(1100, max(300, 28 * n_params)))
    ax_g = Axis(
        fig_g[1, 1];
        xlabel="Year",
        ylabel="Parameter",
        title="Temporal PAWN — $(reef_id) — $(group_names[g])",
        yticks=(1:n_params, string.(sorted_names_g)),
        yreversed=false
    )
    hm_g = heatmap!(
        ax_g, decimal_years, 1:n_params, sorted_pawn_g';
        colormap=:viridis,
        colorrange=(0.0, max(0.1, maximum(sorted_pawn_g)))
    )
    for yr in disturbance_years
        vlines!(ax_g, yr; color=(:red, 0.7), linewidth=1.5, linestyle=:dash)
    end
    Colorbar(fig_g[1, 2], hm_g; label="PAWN mean index", width=14)
    save(
        joinpath(
            ensemble_fig_dir, "$(reef_id)_temporal_pawn_$(group_names[g])_heatmap.png"
        ),
        fig_g; px_per_unit=DPI
    )
end
