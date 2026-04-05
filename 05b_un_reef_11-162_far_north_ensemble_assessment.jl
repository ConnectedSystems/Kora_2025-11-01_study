using BlackBoxOptim
using PairPlots

include("common.jl")

param_bounds = (;
    density=(1.0, 10.0),

    # Group proportions
    prop_tab_acro=(0.001, 0.4),
    prop_cor_acro=(0.001, 0.4),
    prop_cor_non_acro=(0.001, 0.4),
    prop_sm_mass=(0.001, 0.4),
    prop_lrg_mass=(0.001, 0.4),

    # Size dist mean
    size_mean_tab_acro=(1.0, 5.0),
    size_mean_cor_acro=(1.0, 5.0),
    size_mean_cor_non_acro=(1.0, 5.0),
    size_mean_sm_mass=(1.0, 5.0),
    size_mean_lrg_mass=(1.0, 5.0),

    # Size dist stdev
    size_stdev_tab_acro=(0.1, 3.0),
    size_stdev_cor_acro=(0.1, 3.0),
    size_stdev_cor_non_acro=(0.1, 3.0),
    size_stdev_sm_mass=(0.1, 3.0),
    size_stdev_lrg_mass=(0.1, 3.0),

    # Growth scalers
    scalers_tab_acro=(0.25, 1.75),
    scalers_cor_acro=(0.25, 1.75),
    scalers_cor_non_acro=(0.25, 1.75),
    scalers_sm_mass=(0.25, 1.75),
    scalers_lrg_mass=(0.25, 1.75),

    # Recruitment
    recruitment=(0.001, 0.2),
    self_seeding=(0.001, 0.2)
)

reef_id = "11-162"
ensemble_dir = joinpath(OUTPUT_DIR, "ensemble", "offshore_north", "11-162")

ensemble_output = deserialize(joinpath(ensemble_dir, "$(reef_id)_ensemble_output.dat"));
un_reef_ensemble = deserialize(joinpath(ensemble_dir, "$(reef_id)_tracked_candidates.dat"));
ensemble_params = hcat(un_reef_ensemble.candidates...);

identifiability_df = parameter_identifiability_metrics(
    ensemble_params, ENSEMBLE_PARAM_NAMES
)
corr_df = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.4
)
corr_df[!, :Correlation] .= round.(corr_df.Correlation; digits=3)
CSV.write("$(ensemble_dir)/$(reef_id)_parameter_correlations.csv", corr_df)

pawn_sa_results = pawn(ensemble_params', un_reef_ensemble.fitnesses, ENSEMBLE_PARAM_NAMES)

f, ax, sp = heatmap(
    pawn_sa_results[sortperm(pawn_sa_results[PAWNᵢ=At(:median)]), :];
    colormap=:viridis,
    colorrange=(-0.1, maximum(pawn_sa_results))
)
ax.xticklabelrotation[] = π / 2

# Unconstrained sensitivity analysis

# Create paths
ensemble_data_dir = joinpath(OUTPUT_DIR, "sensitivity", "offshore_north", "11-162", "ensemble")
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

ensemble_fig_dir = joinpath(FIG_DIR, "sensitivity", "offshore_north", "11-162", "ensemble")
mkpath(ensemble_fig_dir)

if !isfile(fn_unconstrained_samples)
    # Create objective function
    reef_state = ensemble_output.reef_state
    env_conditions = ensemble_output.env_conditions

    # Load data
    reef_df = CSV.read(
        "data/Reef $(reef_id)_Benthic_line_chart_modelled_2025-12-22.csv", DataFrame
    )
    reef_obs = DataFrame(; SAMPLE_DATE=reef_df.report_year, MEAN_LIVE_CORAL=reef_df.median)

    start_year = Year(Date(reef_obs.SAMPLE_DATE[1])).value
    end_year = Year(Date(reef_obs.SAMPLE_DATE[end])).value - 1  # Exclude last year
    sim_year_range = create_simulation_dates(
        start_year, end_year, calib_settings.start_month
    )
    sim_indices, ref_indices, matched_dates = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
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

# Create data paths
fn_constrained_samples = joinpath(
    ensemble_data_dir, "$(reef_id)_constrained_samples.dat"
)
fn_constrained_fitness = joinpath(
    ensemble_data_dir, "$(reef_id)_constrained_fitness.dat"
)
fn_constrained_pawn = joinpath(ensemble_data_dir, "$(reef_id)_constrained_pawn_results.dat")

if !isfile(fn_constrained_samples)
    # Now constrain to parameter ranges that found good fitness using ensemble values
    # This is to explore the identified range in more depth
    param_bounds = extrema.(eachrow(ensemble_params))
    unif_dists = [Uniform(l, u) for (l, u) in param_bounds]
    n_bounds = length(param_bounds)
    cons_samples = Matrix(
        QMC.sample(n, zeros(n_bounds), ones(n_bounds), sample_method)'
    )

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

# Ensemble correlations
identifiability_df = parameter_identifiability_metrics(
    ensemble_params, ENSEMBLE_PARAM_NAMES
)
corr_threshold = parameter_correlation_analysis(
    ensemble_params, ENSEMBLE_PARAM_NAMES; corr_threshold=0.4
)

corr_params = unique(vcat(corr_threshold.Param1, corr_threshold.Param2))

df = DataFrame(Matrix(ensemble_params'), :auto)
rename!(df, ENSEMBLE_PARAM_NAMES)

# Create pairplot of parameters with absolute Pearson correlation > 0.5
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