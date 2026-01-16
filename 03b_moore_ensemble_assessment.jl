using BlackBoxOptim
using PairPlots

include("common.jl")
include("_16071S_config.jl")

reef_id = "16071S"
ensemble_dir = joinpath(OUTPUT_DIR, "ensemble")

ensemble_output = deserialize(joinpath(ensemble_dir, "$(reef_id)_ensemble_output.dat"));
moore_ensemble = deserialize(joinpath(ensemble_dir, "$(reef_id)_tracked_candidates.dat"));
ensemble_params = hcat(moore_ensemble.candidates...);

# ensemble_params = hcat(calibration_output.tracked_candidates...)
# ensemble_res = calibration_output.ensemble_res
parameter_identifiability_metrics(ensemble_params, ENSEMBLE_PARAM_NAMES)
# parameter_correlation_analysis(ensemble_params, ENSEMBLE_PARAM_NAMES)

corr_df = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.6
)
corr_df[!, :Correlation] .= round.(corr_df.Correlation; digits=3)
CSV.write("$(OUTPUT_DIR)/ensemble/$(reef_id)_parameter_correlations.csv", corr_df)

# SA of ensemble
# pawn_sa_results = pawn(ensemble_params', moore_ensemble.fitnesses, ENSEMBLE_PARAM_NAMES)
# f, ax, sp = heatmap(
#     pawn_sa_results[sortperm(pawn_sa_results[PAWNᵢ=At(:median)]).data, :];
#     colormap=:viridis,
#     colorrange=(-0.1, maximum(pawn_sa_results))
# )
# ax.xticklabelrotation[] = π / 2

# Unconstrained sensitivity analysis

# Create paths
ensemble_data_dir = joinpath(OUTPUT_DIR, "sensitivity", "ensemble")
mkpath(ensemble_data_dir)

fn_unconstrained_samples = joinpath(
    ensemble_data_dir, "$(reef_id)_unconstrained_samples.dat"
)
fn_unconstrained_fitness = joinpath(
    ensemble_data_dir, "$(reef_id)_unconstrained_fitness.dat"
)
fn_unconstrained_pawn = joinpath(
    ensemble_data_dir, "$(reef_id)_unconstrained_pawn_results.dat"
)

ensemble_fig_dir = joinpath(FIG_DIR, "sensitivity", "ensemble")
mkpath(ensemble_fig_dir)

if !isfile(fn_unconstrained_samples)
    # Create objective function
    reef_state = ensemble_output.reef_state
    env_conditions = ensemble_output.env_conditions

    reef_df = CSV.read(
        "$(OUTPUT_DIR)/Moore Reef_Manta Tow_line_chart_modelled_2025-12-28.csv", DataFrame
    )
    reef_obs = DataFrame(;
        SAMPLE_DATE=reef_df.report_year,
        MEAN_LIVE_CORAL=reef_df.mean,
        LOWER=reef_df.lower,
        UPPER=reef_df.upper
    )

    start_year = Year(Date(reef_obs.SAMPLE_DATE[3])).value # start 1994
    end_year = Year(Date(reef_obs.SAMPLE_DATE[end])).value

    sim_year_range = create_simulation_dates(
        start_year, end_year, calib_settings.start_month
    )
    sim_indices, ref_indices, matched_dates = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
    )

    # Exclude specified years
    sim_indices, ref_indices, matched_dates = exclude_years_from_indices(
        sim_indices, ref_indices, matched_dates, reef_config.exclude_years
    )

    model_runner = create_objective_function(
        reef_state, env_conditions, reef_obs.MEAN_LIVE_CORAL,
        sim_indices, ref_indices, 900.0f0,
        89, true
    )

    sample_method = SobolSample(; R=OwenScramble(; base=2, pad=32))

    # Create uniformly distributed samples for uncertain parameters
    n = 8192  # 2^13
    unif_dists = [Uniform(l, u) for (l, u) in param_bounds]
    n_bounds = length(param_bounds)
    unc_samples = Matrix(
        QMC.sample(n, zeros(n_bounds), ones(n_bounds), sample_method)'
    )

    # Create Dirichlet sample for population proportion using Gamma normalization trick
    unc_samples[:, 2:6] = gamma_to_dirichlet(unc_samples[:, 2:6])

    # Scale uniform samples for remaining parameters to uniform distributions using the inverse CDF method
    free_vary_cols = [1, 7:23...]
    unc_samples[:, free_vary_cols] = Matrix(
        Distributions.quantile.(unif_dists[[1, 7:23...]], unc_samples[:, free_vary_cols]')'
    )

    serialize(fn_unconstrained_samples, unc_samples)

    @info "Running unconstrained sample"
    unc_fitness_scores = map(x -> model_runner(collect(x)), eachrow(unc_samples))

    # Save results
    serialize(fn_unconstrained_fitness, unc_fitness_scores)

    unc_pawn_sa_results = pawn(unc_samples, unc_fitness_scores, ENSEMBLE_PARAM_NAMES)
    serialize(fn_unconstrained_pawn, unc_pawn_sa_results)

    f, ax, sp = heatmap(
        unc_pawn_sa_results[sortperm(unc_pawn_sa_results[PAWNᵢ=At(:median)]), :];
        colormap=:viridis,
        colorrange=(-0.1, maximum(unc_pawn_sa_results))
    )
    ax.xticklabelrotation[] = π / 2

    save(joinpath(ensemble_fig_dir, "$(reef_id)_unconstrained_sa.png"), f; px_per_unit=DPI)
else
    unc_samples = deserialize(fn_unconstrained_samples)
    unc_fitness_scores = deserialize(fn_unconstrained_fitness)
    unc_pawn_sa_results = deserialize(fn_unconstrained_pawn)
end

# parameter_identifiability_metrics(unc_samples', ENSEMBLE_PARAM_NAMES)
# parameter_correlation_analysis(unc_samples, ENSEMBLE_PARAM_NAMES)

# Now constrain to parameter ranges that found good fitness using ensemble values
# This is to explore the identified range in more depth

# Create data paths
fn_constrained_samples = joinpath(
    ensemble_data_dir, "$(reef_id)_constrained_samples.dat"
)
fn_constrained_fitness = joinpath(
    ensemble_data_dir, "$(reef_id)_constrained_fitness.dat"
)
fn_constrained_pawn = joinpath(ensemble_data_dir, "$(reef_id)_constrained_pawn_results.dat")

if !isfile(fn_constrained_samples)
    param_bounds = extrema.(eachrow(ensemble_params))
    unif_dists = [Uniform(l, u) for (l, u) in param_bounds]
    n_bounds = length(param_bounds)
    cons_samples = Matrix(
        QMC.sample(n, zeros(n_bounds), ones(n_bounds), sample_method)'
    )

    # Create Dirichlet sample for population proportion using Gamma normalization trick
    cons_samples[:, 2:6] = gamma_to_dirichlet(cons_samples[:, 2:6])

    free_vary_cols = [1, 7:23...]
    cons_samples[:, free_vary_cols] = Matrix(
        Distributions.quantile.(unif_dists[[1, 7:23...]], cons_samples[:, free_vary_cols]')'
    )
    serialize(fn_constrained_samples, cons_samples)

    @info "Running ensemble-constrained sample"
    cons_fitness_scores = map(x -> model_runner(collect(x)), eachrow(cons_samples))

    serialize(fn_constrained_fitness, cons_fitness_scores)

    cons_pawn_sa_results = pawn(cons_samples, cons_fitness_scores, ENSEMBLE_PARAM_NAMES)
    serialize(fn_constrained_pawn, cons_pawn_sa_results)

    f, ax, sp = heatmap(
        cons_pawn_sa_results[sortperm(cons_pawn_sa_results[PAWNᵢ=At(:median)]), :];
        colormap=:viridis,
        colorrange=(-0.1, maximum(cons_pawn_sa_results))
    )
    ax.xticklabelrotation[] = π / 2

    save(joinpath(ensemble_fig_dir, "$(reef_id)_constrained_sa.png"), f; px_per_unit=DPI)
else
    cons_samples = deserialize(fn_constrained_samples)
    cons_fitness_scores = deserialize(fn_constrained_fitness)
    cons_pawn_sa_results = deserialize(fn_constrained_pawn)
end

# Set Makie theme for publication
fontsize_theme = Theme(; fontsize=14)
set_theme!(fontsize_theme)

# Ensemble correlations
parameter_identifiability_metrics(ensemble_params, ENSEMBLE_PARAM_NAMES)
corr_threshold = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.6
)

corr_params = unique(vcat(corr_threshold.Param1, corr_threshold.Param2))

df = DataFrame(Matrix(ensemble_params'), :auto)
rename!(df, ENSEMBLE_PARAM_NAMES)

# Create pairplot of parameters with absolute Pearson correlation > threshold
target_df = df[:, corr_params]
f = pairplot(
    target_df => (
        PairPlots.Series(target_df),
        PairPlots.HexBin(; colormap=Makie.cgrad([:transparent, :black])),
        PairPlots.Scatter(; alpha=0.5),
        PairPlots.Contour(),
        PairPlots.MarginDensity(; bandwidth=0.1),
        PairPlots.MarginQuantileText()
        # PairPlots.MarginQuantileLines()
        # PairPlots.Calculation(
        #     cor;
        #     color=:blue,
        #     position=Makie.Point2f(0.2, 0.1)
        # )
    )
)

autolimits!()
resize_to_layout!(f)
sleep(5)
save("$(ensemble_fig_dir)/$(reef_id)_ensemble_corr_param_pairplot.png", f; px_per_unit=DPI)

cons_pawn_sa_results[
    sortperm(cons_pawn_sa_results[PAWNᵢ=At(:median)]; rev=true), At(:median)
].data

most_influential = collect(
    sortperm(cons_pawn_sa_results[PAWNᵢ=At(:median)]; rev=true)[1:10]
)

# Get factor by names
factor_names = collect(collect(cons_pawn_sa_results.factors[most_influential]))

# Create pairplot of influential parameters
target_df = df[:, factor_names]
f = pairplot(
    target_df => (
        PairPlots.Series(target_df),
        PairPlots.HexBin(; colormap=Makie.cgrad([:transparent, :black])),
        PairPlots.Scatter(; alpha=0.5),
        PairPlots.Contour(),
        PairPlots.MarginDensity(; bandwidth=0.1),
        PairPlots.MarginQuantileText()
        # PairPlots.MarginQuantileLines()
        # PairPlots.Calculation(
        #     cor;
        #     color=:blue,
        #     position=Makie.Point2f(0.2, 0.1)
        # )
    );
    labels=Dict(
        f => rich(string(f); fontsize=18) for f in factor_names
    )
)

autolimits!()
resize_to_layout!(f)
sleep(5)  # Ensure figure generates completely before saving
save(
    "$(ensemble_fig_dir)/$(reef_id)_cons_ensemble_sa_param_pairplot.png", f; px_per_unit=DPI
)
