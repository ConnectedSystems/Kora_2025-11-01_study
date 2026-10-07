using BlackBoxOptim
using PairPlots

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "_16071S_config.jl"))

reef_id = "16071S"
ensemble_dir = joinpath(OUTPUT_DIR, "ensemble", "offshore_north", "moore")

ensemble_output = deserialize(joinpath(ensemble_dir, "$(reef_id)_ensemble_output.dat"));
moore_ensemble = load_result(joinpath(ensemble_dir, "$(reef_id)_tracked_candidates.h5"));
ensemble_params = hcat(moore_ensemble.candidates...);

corr_df = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.5
)
corr_df[!, :Correlation] .= round.(corr_df.Correlation; digits=3)
CSV.write(joinpath(ensemble_dir, "$(reef_id)_parameter_correlations.csv"), corr_df)

# ── Ensemble nearest-neighbour diversity ──────────────────────────────────────
# Flags near-duplicate candidates (repeated optimizer convergence to the same
# basin) that the aggregate `effective_ensemble_diversity` mean can hide.
nn_diversity_df = parameter_nearest_neighbour_diversity(
    ensemble_params, ENSEMBLE_PARAM_NAMES, values(param_bounds); near_dup_threshold=0.05
)
nn_diversity_df[!, :nn_distance] .= round.(nn_diversity_df.nn_distance; digits=4)
nn_diversity_df[!, :rel_nn_distance] .= round.(nn_diversity_df.rel_nn_distance; digits=4)
CSV.write(joinpath(ensemble_dir, "$(reef_id)_nn_diversity.csv"), nn_diversity_df)

n_near_dup = sum(nn_diversity_df.near_duplicate)
if n_near_dup > 0
    @warn "$(n_near_dup)/$(nrow(nn_diversity_df)) ensemble candidates are near-duplicates of another candidate" reef_id
end

# ── Create paths ──────────────────────────────────────────────────────────────
ensemble_data_dir = joinpath(OUTPUT_DIR, "sensitivity", "offshore_north", "moore", "ensemble")
mkpath(ensemble_data_dir)

fn_unconstrained_samples = joinpath(
    ensemble_data_dir, "$(reef_id)_unconstrained_samples.h5"
)
fn_unconstrained_fitness = joinpath(
    ensemble_data_dir, "$(reef_id)_unconstrained_fitness.h5"
)
fn_unconstrained_pawn = joinpath(
    ensemble_data_dir, "$(reef_id)_unconstrained_pawn_results.h5"
)

ensemble_fig_dir = joinpath(FIG_DIR, "sensitivity", "offshore_north", "moore", "ensemble")
mkpath(ensemble_fig_dir)

f_nn = plot_nn_diversity_histogram(
    nn_diversity_df, "Nearest-neighbour diversity: $(reef_id)"
)
save(joinpath(ensemble_fig_dir, "$(reef_id)_nn_diversity_histogram.png"), f_nn; px_per_unit=DPI)

# Loaded unconditionally: sim_year_range is referenced by the temporal PAWN
# section regardless of whether cached sensitivity results exist on disk.
reef_df = CSV.read(
    joinpath(OUTPUT_DIR, "Moore Reef_Manta Tow_line_chart_modelled_2025-12-28.csv"), DataFrame
)
reef_obs = DataFrame(;
    SAMPLE_DATE=reef_df.report_year,
    MEAN_LIVE_CORAL=reef_df.mean,
    LOWER=reef_df.lower,
    UPPER=reef_df.upper
)

start_year = Year(Date(reef_obs.SAMPLE_DATE[3])).value
end_year = Year(Date(reef_obs.SAMPLE_DATE[end])).value

sim_year_range = create_simulation_dates(
    start_year, end_year, calib_settings.start_month
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

benthic_estimate = CSV.read(
    joinpath(OUTPUT_DIR, "ecorrap_benthic", "moore_estimate.csv"), DataFrame
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
    89, calib_settings.use_scalers;
    transform_proportions=false
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

    save_result(fn_unconstrained_samples, unc_samples)

    @info "Running unconstrained sample ($(n_threads) threads)"
    unc_fitness_scores = Vector{Float64}(undef, size(unc_samples, 1))
    Threads.@threads :static for i in axes(unc_samples, 1)
        unc_fitness_scores[i] = runner_pool()(collect(unc_samples[i, :]))
    end

    save_result(fn_unconstrained_fitness, unc_fitness_scores)

    unc_pawn_sa_results = pawn(unc_samples, unc_fitness_scores, ENSEMBLE_PARAM_NAMES)
    save_result(fn_unconstrained_pawn, unc_pawn_sa_results)

    f = plot_pawn_heatmap(
        unc_pawn_sa_results, "Unconstrained SA - $(reef_id)";
        stats=[:mean, :std], xticklabelrotation=π / 2
    )
    save(joinpath(ensemble_fig_dir, "$(reef_id)_unconstrained_sa.png"), f; px_per_unit=DPI)
else
    unc_samples = load_result(fn_unconstrained_samples)
    unc_fitness_scores = load_result(fn_unconstrained_fitness)
    unc_pawn_sa_results = load_result(fn_unconstrained_pawn)
end

# ── Constrained sensitivity analysis ─────────────────────────────────────────
fn_constrained_samples = joinpath(
    ensemble_data_dir, "$(reef_id)_constrained_samples.h5"
)
fn_constrained_fitness = joinpath(
    ensemble_data_dir, "$(reef_id)_constrained_fitness.h5"
)
fn_constrained_pawn = joinpath(ensemble_data_dir, "$(reef_id)_constrained_pawn_results.h5")

if !isfile(fn_constrained_samples)
    # Not `param_bounds`: that must keep the calibration bounds for the Panel B prior below.
    cons_param_bounds = extrema.(eachrow(ensemble_params))
    unif_dists = [Uniform(l, u) for (l, u) in cons_param_bounds]
    n_bounds = length(cons_param_bounds)
    cons_samples = Matrix(
        QMC.sample(n, zeros(n_bounds), ones(n_bounds), sample_method)'
    )
    cons_samples = Matrix(Distributions.quantile.(unif_dists, cons_samples')')

    # Ensemble candidates store realized proportions, so sample each within the ensemble
    # range and renormalize to sum to one (no Gamma transform; the model receives these
    # values as-is via `transform_proportions=false`).
    cons_samples[:, 2:6] ./= sum(cons_samples[:, 2:6]; dims=2)
    save_result(fn_constrained_samples, cons_samples)

    @info "Running ensemble-constrained sample ($(n_threads) threads)"
    cons_fitness_scores = Vector{Float64}(undef, size(cons_samples, 1))
    Threads.@threads :static for i in axes(cons_samples, 1)
        cons_fitness_scores[i] = runner_pool()(
            collect(cons_samples[i, :])
        )
    end

    save_result(fn_constrained_fitness, cons_fitness_scores)

    cons_pawn_sa_results = pawn(cons_samples, cons_fitness_scores, ENSEMBLE_PARAM_NAMES)
    save_result(fn_constrained_pawn, cons_pawn_sa_results)

    f = plot_pawn_heatmap(
        cons_pawn_sa_results, "Constrained SA - $(reef_id)";
        stats=[:mean, :std], xticklabelrotation=π / 2
    )
    save(joinpath(ensemble_fig_dir, "$(reef_id)_constrained_sa.png"), f; px_per_unit=DPI)
else
    cons_samples = load_result(fn_constrained_samples)
    cons_fitness_scores = load_result(fn_constrained_fitness)
    cons_pawn_sa_results = load_result(fn_constrained_pawn)
end

# ── Export PAWN results as CSV ────────────────────────────────────────────────
function _pawn_to_csv(results::AbstractDimArray, path::String)
    factors = string.(collect(dims(results, 1)))
    stats   = collect(dims(results, 2))
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
identifiability_df_full = parameter_identifiability_metrics(ensemble_params, ENSEMBLE_PARAM_NAMES)
identifiability_df_full[!, :CV] .= round.(identifiability_df_full.CV; digits=3)
identifiability_df_full[!, :range_ratio] .= round.(identifiability_df_full.range_ratio; digits=3)
identifiability_df_full[!, :MAD] .= round.(identifiability_df_full.MAD; digits=3)
identifiability_df_full[!, :rMAD] .= round.(identifiability_df_full.rMAD; digits=3)
CSV.write(joinpath(ensemble_data_dir, "$(reef_id)_parameter_identifiability_full.csv"), identifiability_df_full)
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
cons_pawn_sa_results[
    sortperm(cons_pawn_sa_results[PAWNᵢ=At(:mean)]; rev=true), At(:mean)
].data

most_influential = collect(
    sortperm(cons_pawn_sa_results[PAWNᵢ=At(:mean)]; rev=true)[1:10]
)

factor_names = Array(dims(cons_pawn_sa_results, :factors))[most_influential]

target_df = df[:, factor_names]

# Trimmed set (top 5, plus colony density) for the main-text figure; the full
# top-10 pairplot is retained as a supplementary figure (see fig_full below)
# since the complete 10x10 grid is too dense to draw a conclusion from at a
# glance. Density is included for its ecological relevance regardless of PAWN
# rank; it is only appended if not already among the top 5 (its PAWN rank
# varies by reef — e.g. it is ~13th at Moore but can be higher elsewhere).
highlight_factor_names = factor_names[1:5]
if Symbol("Density") ∉ highlight_factor_names
    highlight_factor_names = vcat(highlight_factor_names, [Symbol("Density")])
end
highlight_df = df[:, highlight_factor_names]

# ── Uncertainty reduction: prior vs. posterior ────────────────────────────────
# ensemble_params is (n_params × n_candidates); need (n_candidates × n_params)
posterior_samples = Matrix(ensemble_params')

calib_prior_samples = copy(unc_samples)
latent_proportions = Matrix(QMC.sample(n, zeros(5), ones(5), SobolSample())')
proportion_dists = [Uniform(l, u) for (l, u) in values(param_bounds)[2:6]]
latent_proportions = Matrix(Distributions.quantile.(proportion_dists, latent_proportions')')
calib_prior_samples[:, 2:6] = mapslices(gamma_to_dirichlet, latent_proportions; dims=2)

σ_prior     = vec(std(calib_prior_samples; dims=1))
σ_posterior = vec(std(posterior_samples; dims=1))
rel_reduction = (σ_prior .- σ_posterior) ./ σ_prior

# Top-3 most-sensitive (PAWN) parameters
top3_sensitive = most_influential[1:3]

function plot_uncertainty_panel!(ax, pidx, prior_samples, posterior_samples, rel_reduction)
    unc_v  = prior_samples[:,     pidx]
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
        fontsize=22,
        color=:black
    )
end

# ── Combined figure: constrained SA pairplot (a) + uncertainty reduction (b) ──
# Panel A shows the top-5 most influential factors only. The full top-10
# pairplot (fig_full, below) is provided as a supplementary figure: with 10
# factors the 10x10 grid (45 panels) is too dense to draw a conclusion from
# in the main text.
# Panel A widened/heightened from the original 2000x1950 (5 factors) to fit
# the 6th factor (Density) added alongside the top-5 without crowding labels.
fig_combined = Figure(; size=(2400, 3150))

pp_grid = pairplot(
    fig_combined[1, 1],
    highlight_df => (
        PairPlots.Series(highlight_df),
        PairPlots.HexBin(; colormap=Makie.cgrad([:transparent, :black])),
        PairPlots.Scatter(; alpha=0.5),
        PairPlots.Contour(),
        PairPlots.MarginDensity(; bandwidth=0.1),
        PairPlots.MarginQuantileText(; fontsize=32)
    );
    labels=Dict(fn => rich(string(fn); fontsize=28) for fn in highlight_factor_names),
    bodyaxis=(; xticklabelsize=20, yticklabelsize=20),
    diagaxis=(; xticklabelsize=20, yticklabelsize=20)
)

# Hexbin shading/contours are scaled per-panel (each panel's own min-to-max
# density), so a single absolute colorbar would misrepresent panels against
# each other. The colorbar below anchors the relative low-to-high convention
# used throughout instead of a shared count scale.
# Placed inside the empty upper-right triangle of the (lower-triangular)
# pairplot grid, rather than as a separate tall/thin column alongside it.
n_highlight = length(highlight_factor_names)
Colorbar(
    pp_grid[1:2, (n_highlight - 3):(n_highlight - 1)];
    vertical=false,
    valign=:top,
    colormap=Makie.cgrad([:transparent, :black]),
    limits=(0, 1),
    ticks=([0, 1], ["Low", "High"]),
    label="Relative ensemble density (per panel)",
    height=30,
    labelsize=28,
    ticklabelsize=24,
)

gl_unc = fig_combined[2, 1] = GridLayout()

# "Most influential factors (PAWN)" describes what the three panels below are
# (the top-3 PAWN-ranked parameters), not a shared y-axis unit (the y-axis is
# "Density") — so it is set as a title above the row rather than a rotated
# axis-style label to the side.
Label(gl_unc[0, 1:3]; text="Most influential factors (PAWN)", fontsize=26, font=:bold, tellwidth=false)

for (col, pidx) in enumerate(top3_sensitive)
    ax = Axis(
        gl_unc[2, col];
        xlabel=string(ENSEMBLE_PARAM_NAMES[pidx]),
        ylabel=col == 1 ? "Density" : "",
        xlabelsize=28,
        ylabelsize=28,
        xticklabelsize=24,
        yticklabelsize=24,
    )
    plot_uncertainty_panel!(ax, pidx, calib_prior_samples, posterior_samples, rel_reduction)
end

Legend(gl_unc[1, 1:3],
    [PolyElement(; color=(:steelblue, 0.6)), PolyElement(; color=(:orangered, 0.6))],
    ["Prior (calibration bounds)", "Posterior (calibrated ensemble)"];
    orientation=:horizontal,
    framevisible=false,
    labelsize=24,
)

Label(fig_combined[1, 1, TopLeft()], "(A)"; fontsize=32, font=:bold, padding=(4, 0, 4, 0))
Label(fig_combined[2, 1, TopLeft()], "(B)"; fontsize=32, font=:bold, padding=(4, 0, 4, 0))

rowsize!(fig_combined.layout, 1, 2340)
rowsize!(fig_combined.layout, 2, 350)

sleep(5)
save(
    joinpath(ensemble_fig_dir, "$(reef_id)_cons_ensemble_sa_param_pairplot.png"),
    fig_combined; px_per_unit=DPI
)

# ── Supplementary figure: full top-10 pairplot (Panel A only) ────────────────
# Sized up from the original 2000x2000 to give the larger axis-label fontsize
# (below) room; at 2000x2000 the 10-column grid's bottom/left labels overlapped.
fig_full = Figure(; size=(2600, 2600))

pp_grid_full = pairplot(
    fig_full[1, 1],
    target_df => (
        PairPlots.Series(target_df),
        PairPlots.HexBin(; colormap=Makie.cgrad([:transparent, :black])),
        PairPlots.Scatter(; alpha=0.5),
        PairPlots.Contour(),
        PairPlots.MarginDensity(; bandwidth=0.1),
        PairPlots.MarginQuantileText(; fontsize=26)
    );
    labels=Dict(fn => rich(string(fn); fontsize=22) for fn in factor_names),
    bodyaxis=(; xticklabelsize=16, yticklabelsize=16),
    diagaxis=(; xticklabelsize=16, yticklabelsize=16)
)

# Placed inside the empty upper-right triangle of the pairplot grid, rather
# than as a separate tall/thin column alongside it.
n_full = length(factor_names)
Colorbar(
    pp_grid_full[1:3, (n_full - 4):(n_full - 2)];
    vertical=false,
    valign=:top,
    colormap=Makie.cgrad([:transparent, :black]),
    limits=(0, 1),
    ticks=([0, 1], ["Low", "High"]),
    label="Relative ensemble density (per panel)",
    height=30,
    labelsize=26,
    ticklabelsize=22,
)

sleep(5)
save(
    joinpath(ensemble_fig_dir, "$(reef_id)_cons_ensemble_sa_param_pairplot_full.png"),
    fig_full; px_per_unit=DPI
)

# ═══════════════════════════════════════════════════════════════════════════════
# ADDITION 2 — Temporal sensitivity of coral cover to recruitment
#
# The constrained sample (parameters drawn independently within the ensemble ranges)
# is run twice: with recruitment as sampled, and with both recruitment pathways
# (external recruitment and self-seeding) switched off. The two runs of each sample
# share a random stream, so their difference is the cover attributable to
# recruitment. Dummy-normalized PAWN is computed per year on cover with recruitment,
# without it, and on the difference.
#
# The calibrated ensemble is not used here: its two recruitment parameters trade off
# against each other (strong negative correlation), and 250 members leave too few
# samples per PAWN slice, so neither pathway's influence can be resolved from it.
#
# Robustness checks (CSV only): a 50-dummy 95% null as a stricter bar than the single
# built-in dummy, and PAWN repeated on the samples meeting the calibration fitness
# threshold.
# ═══════════════════════════════════════════════════════════════════════════════

recruit_years = collect(1994:2023)
recruit_seed = 89
recruit_scenarios = [:with, :without, :difference]
recruit_rows = findall(in(["External Recruitment", "Self-seeding"]), ENSEMBLE_PARAM_NAMES)
disturbance_years = Float64.(reef_config.disturbance_years)
disturbance_label = "2017 Bleaching event"

fn_recruit_ts_with = joinpath(ensemble_data_dir, "$(reef_id)_recruitment_ts_with.h5")
fn_recruit_ts_without = joinpath(ensemble_data_dir, "$(reef_id)_recruitment_ts_without.h5")
fn_recruit_pawn = joinpath(ensemble_data_dir, "$(reef_id)_recruitment_temporal_pawn.h5")

# ── Paired model runs ─────────────────────────────────────────────────────────
if !isfile(fn_recruit_ts_with) || !isfile(fn_recruit_ts_without)
    @info "Running constrained sample with and without recruitment " *
        "($(size(cons_samples, 1)) × 2 runs, $(Threads.nthreads()) threads)"
    recruit_ts_runner = create_timeseries_function(
        reef_state, env_conditions, reef_config.area, recruit_seed, true;
        transform_proportions=false
    )
    ts_with, ts_without = recruitment_paired_runs(recruit_ts_runner, cons_samples, recruit_seed)
    save_result(fn_recruit_ts_with, ts_with)
    save_result(fn_recruit_ts_without, ts_without)
else
    ts_with = load_result(fn_recruit_ts_with)
    ts_without = load_result(fn_recruit_ts_without)
end
@assert size(ts_with, 2) == length(recruit_years) "recruit_years does not match the simulation length"

recruit_ts = Dict(
    :with => ts_with, :without => ts_without, :difference => ts_with .- ts_without
)

# ── Temporal PAWN per scenario ────────────────────────────────────────────────
if !isfile(fn_recruit_pawn)
    @info "Computing temporal PAWN with and without recruitment"
    Random.seed!(42)  # reproducible PAWN dummy factor
    recruit_pawn = DimArray(
        cat(
            [parent(temporal_pawn(cons_samples, recruit_ts[s], ENSEMBLE_PARAM_NAMES, recruit_years))
             for s in recruit_scenarios]...;
            dims=3
        ),
        (
            Dim{:factors}(Symbol.(ENSEMBLE_PARAM_NAMES)),
            Dim{:year}(recruit_years),
            Dim{:scenario}(recruit_scenarios)
        )
    )
    save_result(fn_recruit_pawn, recruit_pawn)
else
    recruit_pawn = load_result(fn_recruit_pawn)
end

# Long-format table (one row per scenario × year × parameter) for reporting
recruit_pawn_df = DataFrame(;
    scenario=String[], year=Int[], parameter=String[], pawn=Float64[], rank=Int[]
)
for s in recruit_scenarios, yr in recruit_years
    vals = recruit_pawn[scenario=At(s), year=At(yr)].data
    ranks = invperm(sortperm(vals; rev=true))
    for (j, pname) in enumerate(ENSEMBLE_PARAM_NAMES)
        push!(recruit_pawn_df, (string(s), yr, pname, round(vals[j]; digits=4), ranks[j]))
    end
end
CSV.write(
    joinpath(ensemble_data_dir, "$(reef_id)_recruitment_temporal_pawn.csv"), recruit_pawn_df
)

# ── Robustness checks ─────────────────────────────────────────────────────────
# (1) Factors above the 95th percentile of 50 random dummies (raw mean-KS PAWN), and
#     whether each recruitment pathway clears that bar.
# (2) Recruitment ranks when PAWN is restricted to the best-fitting 10% of samples. No
#     constrained sample reaches the calibration fitness threshold (good fits need the
#     correlated parameter combinations the ensemble found), so the best decile is the
#     closest available subset. Parameters are correlated within it, so it is a check
#     on direction, not a replacement analysis.
fit_cutoff = quantile(cons_fitness_scores, 0.1)
fitness_ok = cons_fitness_scores .<= fit_cutoff
n_fit = count(fitness_ok)
@info "Best-fitting decile of the constrained sample: $(n_fit) samples, fitness <= " *
    "$(round(fit_cutoff; digits=3)) (calibration threshold: $(opt_config.fitness_threshold))"

robust_df = DataFrame(;
    scenario=String[], year=Int[], null95=Float64[], n_above_null95=Int[],
    ext_rank=Int[], ext_above_null95=Bool[], self_rank=Int[], self_above_null95=Bool[],
    n_best_decile=Int[], ext_rank_best_decile=Union{Int,Missing}[], self_rank_best_decile=Union{Int,Missing}[]
)
Random.seed!(42)
for s in recruit_scenarios, (t, yr) in enumerate(recruit_years)
    y = Float64.(recruit_ts[s][:, t])
    std(y) == 0 && continue  # e.g. the difference in the first year

    thr, raw = pawn_dummy_null(cons_samples, y)
    ranks = invperm(sortperm(raw; rev=true))

    ranks_fit = (missing, missing)
    y_fit = y[fitness_ok]
    if std(y_fit) > 0
        Si_fit = (
            pawn(cons_samples[fitness_ok, :], y_fit, ENSEMBLE_PARAM_NAMES)[PAWNᵢ=At(:mean)].data
        )
        r_fit = invperm(sortperm(Si_fit; rev=true))
        ranks_fit = (r_fit[recruit_rows[1]], r_fit[recruit_rows[2]])
    end

    push!(robust_df, (
        string(s), yr, round(thr; digits=4), count(>(thr), raw),
        ranks[recruit_rows[1]], raw[recruit_rows[1]] > thr,
        ranks[recruit_rows[2]], raw[recruit_rows[2]] > thr,
        n_fit, ranks_fit...
    ))
end
CSV.write(
    joinpath(ensemble_data_dir, "$(reef_id)_recruitment_temporal_robustness.csv"), robust_df
)

# ── Figures ───────────────────────────────────────────────────────────────────
fig_title = "$(reef_config.reef_name): constrained sample (n = $(size(cons_samples, 1)))"
save(
    joinpath(ensemble_fig_dir, "$(reef_id)_recruitment_trajectories.png"),
    plot_recruitment_trajectories(
        recruit_years, ts_with, ts_without; disturbance_years=disturbance_years,
        disturbance_label=disturbance_label, title=fig_title
    );
    px_per_unit=DPI
)
save(
    joinpath(ensemble_fig_dir, "$(reef_id)_recruitment_temporal_pawn_heatmaps.png"),
    plot_recruitment_temporal_pawn(
        recruit_pawn; disturbance_years=disturbance_years,
        disturbance_label=disturbance_label, title=fig_title
    );
    px_per_unit=DPI
)

# Panel C alone (recruitment-attributable cover), for presentations
save(
    joinpath(ensemble_fig_dir, "$(reef_id)_recruitment_temporal_pawn_difference_heatmap.png"),
    plot_recruitment_temporal_pawn(
        recruit_pawn; scenarios=[:difference], disturbance_years=disturbance_years,
        disturbance_label=disturbance_label,
        title=fig_title
    );
    px_per_unit=DPI
)
