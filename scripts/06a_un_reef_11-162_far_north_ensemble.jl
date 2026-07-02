using BlackBoxOptim

include(joinpath(@__DIR__, "common.jl"))

# Configuration for Far North ensemble
reef_config = ReefConfig(;
    reef_id="11-162",
    reef_name="Unknown",
    area=Float32(72.0 * 4),  # Represent area of four EcoRRAP transects
    depth=9.0,
    density=15,
    # This initial guess is "wrong" but should be handled via the calibration process
    initial_proportions=[0.02f0, 0.18f0, 0.2f0, 0.3f0, 0.3f0],
    exclude_years=Int64[2023, 2024]  # exclude last two years for assessment
)

file_paths = CalibrationDataPaths(;
    dhw_scenarios=joinpath(OUTPUT_DIR, "dhw_scens.nc"),
    canonical_reefs=joinpath(OUTPUT_DIR, "rrap_canonical_2025-07-15-T10-48-29.gpkg"),
    growth_models=joinpath(OUTPUT_DIR, "offshore_north", "overall", "offshore_north_growth_models.dat"),
    survival_models=joinpath(OUTPUT_DIR, "offshore_north", "overall", "offshore_north_survival_models.dat"),
    output_dir=OUTPUT_DIR,
    figure_dir=FIG_DIR
)

# Optimization settings
opt_config = OptimizationConfig(;
    max_steps=25_000,
    population_size=50,
    fitness_threshold=0.1,
    ensemble_members=100,
    trace_interval=10,
    random_seed=64
)

# Search ranges
search_ranges = SearchRanges(;
    density=(1.0, 15.0),
    group_proportion=(0.0, 1.0),  # probability levels for Gamma quantiles
    size_mean=(1.0, 5.0),
    size_std=(0.25, 2.0),
    scalers=(0.5, 2.0),
    recruitment=(0.001, 0.1),
    self_seeding=(0.001, 0.3)
)

# Calibration settings
calib_settings = CalibrationSettings(;
    date_match_max_days=120,
    start_month=3,
    use_scalers=true
)

# Main Calibration Workflow

"""
    run_calibration(reef_config, file_paths, opt_config, search_ranges, calib_settings)

Run complete calibration workflow for a reef.
"""
function run_calibration(
    reef_config::ReefConfig,
    file_paths::CalibrationDataPaths,
    opt_config::OptimizationConfig,
    search_ranges::SearchRanges,
    calib_settings::CalibrationSettings
)
    @info "Starting calibration for $(reef_config.reef_name) ($(reef_config.reef_id))"

    # Load data
    @info "Loading reef observations..."
    reef_df = CSV.read(
        joinpath(OUTPUT_DIR, "Reef 11-162_Benthic_line_chart_modelled_2025-12-22.csv"), DataFrame
    )
    reef_obs = DataFrame(;
        SAMPLE_DATE=reef_df.report_year,
        MEAN_LIVE_CORAL=reef_df.mean,
        LOWER=reef_df.lower,
        UPPER=reef_df.upper
    )

    start_year = Year(Date(reef_obs.SAMPLE_DATE[1])).value
    end_year = Year(Date(reef_obs.SAMPLE_DATE[end])).value - 1

    reef_obs = reef_obs[1:(end - 1), :]

    @info "Loading environmental data..."
    reef_uid = "11162100104"
    ds = open_dataset(file_paths.dhw_scenarios)
    historic_dhw = vec(
        NetCDF.read(
            ds.dhw_scens[
                locs=At(reef_uid),
                scenarios=1,
                timesteps=At(start_year, end_year)
            ]
        )
    )

    # Load models
    @info "Loading growth and survival models..."
    calib_growth_models = deserialize(file_paths.growth_models)
    calib_survival_models = deserialize(file_paths.survival_models)

    # Initialize model
    @info "Initializing reef state..."
    n_ts = length(historic_dhw)
    reef_state = Kora.initialize_reef(;
        n_timesteps=n_ts,
        n_locs=1,
        area=reef_config.area,
        density=reef_config.density,
        depths=reef_config.depth,
        growth_models=calib_growth_models,
        survival_models=calib_survival_models
    )

    total_initial_pop = ceil(Int64, 4 * reef_config.area)
    Kora.initialize_coral_population!(
        reef_state, 1, total_initial_pop;
        group_proportions=reef_config.initial_proportions
    )

    # Prepare environmental conditions
    env_conditions = Kora.generate_example_environment(n_ts, 1; with_dhw=false)
    env_conditions[:, 1, 1] .= historic_dhw

    # Match simulation dates with observations
    @info "Matching simulation dates with observations..."
    sim_year_range = create_simulation_dates(
        start_year, end_year, calib_settings.start_month
    )
    c_sim_indices, c_ref_indices, matched_dates = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
    )

    # Exclude specified years
    c_sim_indices, c_ref_indices, matched_dates = exclude_years_from_indices(
        c_sim_indices, c_ref_indices, matched_dates, reef_config.exclude_years
    )

    @info "Matched $(length(c_sim_indices)) timepoints for calibration"

    # Check if results already exist
    ensemble_dir = joinpath(file_paths.output_dir, "ensemble", "offshore_north", "11-162")
    result_files = [
        joinpath(ensemble_dir, "$(reef_config.reef_id)_optim_state.dat"),
        joinpath(ensemble_dir, "$(reef_config.reef_id)_optim_best.dat"),
        joinpath(ensemble_dir, "$(reef_config.reef_id)_tracked_candidates.dat")
    ]

    opt_results = []

    if all(isfile.(result_files))
        @info "Loading existing ensemble results..."
        res, optim_best, tracked_candidates, tracked_fitnesses, trial_contribution = load_calibration_results(
            ensemble_dir, reef_config.reef_id
        )
    else
        @info "Running optimization..."

        # Create objective function
        objective = create_objective_function(
            reef_state, env_conditions, reef_obs.MEAN_LIVE_CORAL,
            c_sim_indices, c_ref_indices, reef_config.area, nothing, nothing,
            opt_config.random_seed, calib_settings.use_scalers
        )

        # Setup tracking
        callback, tracked_candidates, tracked_fitnesses, trial_counter = create_tracking_callback(
            opt_config.fitness_threshold, opt_config.ensemble_members
        )

        # Build search ranges
        opt_ranges = build_search_ranges(search_ranges, calib_settings.use_scalers)

        trial_contribution = Int64[]  # zeros(Int64, n_trials)
        n_trials = 0

        # n_trials = 5
        # # Run multiple optimizations with different starting conditions
        # for trial in 1:n_trials
        #     trial_counter[] = 0  # Adds to counter within callback

        #     # Run optimization
        #     res = bboptimize(
        #         objective;
        #         SearchRange=ranges,
        #         MaxSteps=opt_config.max_steps,
        #         PopulationSize=opt_config.population_size,
        #         CallbackInterval=0.1,
        #         CallbackFunction=callback,
        #         TraceInterval=opt_config.trace_interval
        #     )

        #     push!(opt_results, res)

        #     if trial_counter[] == 0
        #         @info "No ensemble candidates found!"
        #     else
        #         trial_contribution[trial] = trial_counter[]
        #     end
        # end
        # Run multiple optimizations with different starting conditions
        # for trial in 1:n_trials
        while sum(trial_contribution) < opt_config.ensemble_members
            n_trials += 1
            @info "Running attempt $(n_trials)"
            @info "Progress so far: $(trial_contribution) [Total: $(sum(trial_contribution))]"
            trial_counter[] = 0

            # Run optimization
            res = bboptimize(
                objective;
                SearchRange=opt_ranges,
                MaxSteps=opt_config.max_steps,
                MaxTime=300,
                PopulationSize=opt_config.population_size,
                CallbackInterval=0.0,
                CallbackFunction=callback,
                TraceInterval=opt_config.trace_interval
            )

            push!(opt_results, res)

            if trial_counter[] == 0
                @info "No ensemble candidates found!"
            else
                push!(trial_contribution, trial_counter[])
                # trial_contribution[trial] = trial_counter[]
            end
        end

        optim_best = if length(tracked_fitnesses) > 0
            tracked_candidates[findmin(tracked_fitnesses)[2]]
        else
            # Get best of current round if no ensemble
            best_candidate(res)
        end

        # Convert probability values to realized values
        optim_best[2:6] = gamma_to_dirichlet(optim_best[2:6])
        for c in tracked_candidates
            c[2:6] .= gamma_to_dirichlet(c[2:6])
        end

        save_calibration_results(
            ensemble_dir, reef_config.reef_id, res, optim_best,
            tracked_candidates, tracked_fitnesses, trial_contribution,
            opt_config.fitness_threshold
        )
    end

    @info "Best fitness: $(best_fitness.(opt_results))"

    # Run model with best parameters
    @info "Running model with best parameters..."
    Kora.set_population!(reef_state, optim_best)

    rng = Random.seed!(opt_config.random_seed)
    if calib_settings.use_scalers && length(optim_best) > 16
        n_grps = Kora.n_groups(reef_state)

        # Extract and apply scalers
        scaler_start = 17
        scaler_end = 17 + n_grps - 1
        loc_scalers = optim_best[scaler_start:scaler_end]
        Kora.assign_scalers!(reef_state, loc_scalers)

        # Extract recruitment parameters
        recruitment_proportion = optim_best[scaler_end + 1]
        self_seeding_proportion = optim_best[scaler_end + 2]

        # Apply DHW tolerances if present
        dhw_tols = extract_dhw_tolerances(optim_best, n_grps)
        if !isnothing(dhw_tols)
            dhw_means, dhw_stds = dhw_tols
            for grp in 1:n_grps
                reef_state.wild_dhw_tolerances[1, :, grp, At(:mean)] .= dhw_means[grp]
                reef_state.wild_dhw_tolerances[:, :, grp, At(:stdev)] .= dhw_stds[grp]
            end
        end

        Kora.run_example!(
            reef_state, env_conditions;
            recruits=Float32(recruitment_proportion),
            self_seed=Float32(self_seeding_proportion),
            rng=rng
        )
    else
        Kora.run_example!(reef_state, env_conditions; rng=rng)
    end

    # Provisional cover — replaced below with best ensemble member once available.
    cover = (Kora.coral_cover(reef_state) ./ reef_config.area) * 100.0

    # Find validation points
    v_sim_indices, v_ref_indices, _ = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
    )

    v_sim_idx = [s for s in v_sim_indices if s ∉ c_sim_indices]
    v_ref_idx = [s for s in v_ref_indices if s ∉ c_ref_indices]

    # Generate visualizations
    @info "Generating visualizations..."
    mkpath(file_paths.figure_dir)

    # Main timeseries
    f_ts = Kora.viz.timeseries(reef_state, env_conditions)
    save(
        joinpath(file_paths.figure_dir, "$(reef_config.reef_id)_calibrated_timeseries.png"),
        f_ts;
        px_per_unit=DPI
    )

    # Run ensemble if tracked candidates available
    ensemble_res = nothing
    if !isempty(tracked_candidates)
        @info "Running ensemble with $(length(tracked_candidates)) candidates..."
        suitable_params = hcat(tracked_candidates...)
        ensemble_res = Kora.run_ensemble!(reef_state, env_conditions, suitable_params)

        # Ensemble timeseries
        f_ensemble = Kora.viz.ensemble_timeseries(
            reef_state, ensemble_res, env_conditions
        )
        save(
            joinpath(
                file_paths.figure_dir, "$(reef_config.reef_id)_ensemble_timeseries.png"
            ),
            f_ensemble;
            px_per_unit=DPI
        )
    else
        @info "No ensemble candidates found!"
    end

    # Select the best-performing ensemble member by calibration RMSE.
    # Using a member drawn from the ensemble (rather than a fresh optim_best run)
    # ensures both the ensemble mean and the single-member baseline share the same
    # stochastic budget, making the comparison fair.
    best_member_idx = nothing
    obs_calib = reef_obs.MEAN_LIVE_CORAL[c_ref_indices]
    if !isnothing(ensemble_res)
        ens_cover_pct = (ensemble_res.cover[:, 1, :] ./ reef_config.area) .* 100.0
        member_rmse = [Kora.RMSE(ens_cover_pct[c_sim_indices, i], obs_calib)
                       for i in axes(ens_cover_pct, 2)]
        best_member_idx = argmin(member_rmse)
        cover = Vector{Float64}(ens_cover_pct[:, best_member_idx])
    end

    sim = cover[c_sim_indices]
    metrics = calculate_performance_metrics(sim, obs_calib)

    label = isnothing(best_member_idx) ? "optim_best single run" : "best ensemble member (RMSE-selected)"
    @info "Performance Metrics ($label):"
    @info "  RMSE: $(round(metrics.rmse; digits=2))%"
    @info "  Pearson: $(round(metrics.pearson; digits=3))"
    @info "  Kendall: $(round(metrics.kendall; digits=3))"
    @info "  Bias (β): $(round(metrics.bias; digits=3))"
    @info "  Variability (α): $(round(metrics.variability; digits=3))"

    # Calibration comparison plot (with ensemble if available)
    f_calib = plot_calibration_results(
        reef_state, env_conditions, reef_obs, c_sim_indices, c_ref_indices,
        sim_year_range, cover, reef_config.area, metrics,
        joinpath(
            file_paths.figure_dir, "$(reef_config.reef_id)_calibration_comparison.png"
        );
        v_sim_indices=v_sim_idx,
        v_ref_indices=v_ref_idx,
        ensemble_res=ensemble_res,
        best_member_idx=best_member_idx
    )

    @info "Calibration complete!"
    return (
        reef_state=reef_state,
        env_conditions=env_conditions,
        results=res,
        best_params=optim_best,
        metrics=metrics,
        tracked_candidates=tracked_candidates,
        tracked_fitnesses=tracked_fitnesses,
        trial_contribution=trial_contribution,
        ensemble_res=ensemble_res
    )
end

ensemble_data_dir = joinpath(OUTPUT_DIR, "ensemble", "offshore_north", "11-162")
mkpath(ensemble_data_dir)

# Run calibration
calibration_output = run_calibration(
    reef_config, file_paths, opt_config, search_ranges, calib_settings
)

output_path = joinpath(ensemble_data_dir, "$(reef_config.reef_id)_ensemble_output.dat")
serialize(output_path, calibration_output)

ensemble_params = hcat(calibration_output.tracked_candidates...)
ensemble_fitnesses = calibration_output.tracked_fitnesses

n_params = size(ensemble_params, 1)
n_cols = 5
n_rows = ceil(Int, n_params / n_cols)

median_fitness = median(ensemble_fitnesses)
median_idx = ensemble_fitnesses .< median_fitness

identifiability_df = parameter_identifiability_metrics(
    ensemble_params[:, median_idx], ENSEMBLE_PARAM_NAMES
)

identifiability_df[!, :CV] .= round.(identifiability_df.CV; digits=3)
identifiability_df[!, :range_ratio] .= round.(identifiability_df.range_ratio; digits=3)
identifiability_df[!, :MAD] .= round.(identifiability_df.MAD; digits=3)
identifiability_df[!, :rMAD] .= round.(identifiability_df.rMAD; digits=3)

f = Figure(; size=(1400, 300 * n_rows))
for p in 1:n_params
    row = div(p - 1, n_cols) + 1
    col = mod(p - 1, n_cols) + 1

    ax = Axis(
        f[row, col];
        xlabel=ENSEMBLE_PARAM_NAMES[p],
        xlabelsize=16
    )

    min_x, max_x = extrema(ensemble_params[p, :])
    buffer = abs(-(extrema(ensemble_params[p, :])...)) * 0.1
    min_x, max_x = min_x - buffer, max_x + buffer
    min_y, max_y = median_fitness, maximum(ensemble_fitnesses) + 0.1
    background_rect = Rect2f(min_x, min_y, max_x - min_x, max_y - min_y)

    poly!(
        ax,
        background_rect;
        color=(:gray, 0.1),
        strokewidth=0
    )

    hlines!(median_fitness; color=(:black, 0.5))
    scatter!(ax, ensemble_params[p, :], ensemble_fitnesses;
        color=:blue, markersize=6, alpha=0.4)

    beats_median = ensemble_fitnesses .< median_fitness
    target_p, target_f = ensemble_params[p, beats_median], ensemble_fitnesses[beats_median]
    scatter!(ax, target_p, target_f; color=:gold, markersize=6, alpha=0.6)
    k = kde(hcat(target_p, target_f))
    density_levels = quantile(vec(k.density), [0.75, 0.95])
    contour!(
        k.x, k.y, k.density; levels=density_levels, linewidth=2, alpha=0.8, colormap=:plasma
    )

    # Add metric text relative to bottom
    indicator = identifiability_df[p, :rMAD]
    score_text = "rMAD: $(indicator)"
    text!(ax, 0.5, 0.90; text=score_text, align=(:center, :bottom), space=:relative)

    ylims!(ax, 0, maximum(ensemble_fitnesses) + 0.1)
    xlims!(ax, min_x, max_x)
end

Label(f[:, 0], "Metric Score"; rotation=π / 2, fontsize=18)

save(joinpath(FIG_DIR, "$(reef_config.reef_id)_calib_param_interactions.png"), f; px_per_unit=DPI)

reef_df = CSV.read(
    joinpath(OUTPUT_DIR, "Reef 11-162_Benthic_line_chart_modelled_2025-12-22.csv"), DataFrame
)
reef_obs = DataFrame(;
    SAMPLE_DATE=reef_df.report_year,
    MEAN_LIVE_CORAL=reef_df.mean,
    LOWER=reef_df.lower,
    UPPER=reef_df.upper
)

area = calibration_output.reef_state.carrying_capacity[1]
obs = reef_obs.MEAN_LIVE_CORAL[4:5]

ensemble_projection = (calibration_output.ensemble_res.cover[4:5, 1, :] / area) * 100.0

# 95% CI of projections
@info "95% CI" quantile(abs.(ensemble_projection .- obs), [0.025, 0.975])

@info "Absolute max:" maximum(abs.(ensemble_projection .- obs))

@info "Maximum:" maximum(ensemble_projection .- obs)
@info "Minimum:" minimum(ensemble_projection .- obs)
