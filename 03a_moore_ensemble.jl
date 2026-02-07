using BlackBoxOptim

include("common.jl")
include("_16071S_config.jl")

"""
    run_calibration(reef_config, file_paths, opt_config, search_ranges, calib_settings)

Run complete ensemble workflow for a reef.
"""
function run_calibration(
    reef_config::ReefConfig,
    file_paths::CalibrationDataPaths,
    opt_config::OptimizationConfig,
    search_ranges::NamedTuple, # SearchRanges,
    calib_settings::CalibrationSettings
)
    @info "Starting ensemble search for $(reef_config.reef_name) ($(reef_config.reef_id))"

    # Load data
    @info "Loading reef observations..."
    reef_df = CSV.read(
        "$(OUTPUT_DIR)/Moore Reef_Manta Tow_line_chart_modelled_2025-12-28.csv", DataFrame
    )
    reef_obs = DataFrame(;
        SAMPLE_DATE=reef_df.report_year,
        MEAN_LIVE_CORAL=reef_df.mean,
        LOWER=reef_df.lower,
        UPPER=reef_df.upper
    )

    @info "Loading environmental data..."
    reef_uid = get_reef_uid(file_paths.canonical_reefs, reef_config.reef_id)

    start_year = Year(Date(reef_obs.SAMPLE_DATE[3])).value # start 1994
    end_year = Year(Date(reef_obs.SAMPLE_DATE[end])).value

    historic_dhw = load_historical_dhw(
        file_paths.dhw_scenarios, reef_uid, start_year, end_year
    )

    # Load models
    @info "Loading growth and survival models..."
    calib_growth_models = deserialize(file_paths.growth_models)
    calib_survival_models = deserialize(file_paths.survival_models)

    # Initialize model
    @info "Initializing reef state..."
    n_ts = length(historic_dhw)
    reef_state = CoralFlow.initialize_reef(;
        n_timesteps=n_ts,
        n_locs=1,
        area=reef_config.area,
        density=reef_config.density,
        depths=reef_config.depth,
        growth_models=calib_growth_models,
        survival_models=calib_survival_models
    )

    # This is for initialization - the evaluated configuration will be determined
    # by the optimization process
    total_initial_pop = ceil(Int64, 4 * reef_config.area)
    CoralFlow.initialize_coral_population!(
        reef_state, 1, total_initial_pop;
        group_proportions=reef_config.initial_proportions
    )

    # Prepare environmental conditions
    env_conditions = CoralFlow.generate_example_environment(n_ts, 1; with_dhw=false)
    env_conditions[:, 1, 1] .= historic_dhw

    # Match simulation dates with observations
    @info "Matching simulation dates with observations..."
    sim_year_range = create_simulation_dates(
        start_year, end_year, calib_settings.start_month
    )
    c_sim_indices, c_ref_indices, c_matched_dates = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
    )

    # Exclude specified years
    c_sim_indices, c_ref_indices, c_matched_dates = exclude_years_from_indices(
        c_sim_indices, c_ref_indices, c_matched_dates, reef_config.exclude_years
    )

    @info "Matched $(length(c_sim_indices)) timepoints for ensemble search"

    # Check if results already exist
    ensemble_dir = joinpath(file_paths.output_dir, "ensemble")
    result_files = [
        joinpath(ensemble_dir, "$(reef_config.reef_id)_optim_state.dat"),
        joinpath(ensemble_dir, "$(reef_config.reef_id)_optim_best.dat"),
        joinpath(ensemble_dir, "$(reef_config.reef_id)_tracked_candidates.dat")
    ]

    benthic_estimate = CSV.read("data/ecorrap_benthic/moore_estimate.csv", DataFrame)

    year_span = year.(c_matched_dates)
    sim_benthic_years = [year_span .∈ Ref(benthic_estimate.year)][1]
    aligned_years = year_span[sim_benthic_years]
    benthic_data = benthic_estimate[benthic_estimate.year .∈ aligned_years, :]

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
            c_sim_indices, c_ref_indices, reef_config.area, sim_benthic_years, benthic_data,
            opt_config.random_seed, calib_settings.use_scalers
        )

        # Setup tracking
        callback, tracked_candidates, tracked_fitnesses, trial_counter = create_tracking_callback(
            opt_config.fitness_threshold, opt_config.ensemble_members
        )

        # Build search ranges
        # opt_ranges = build_search_ranges(search_ranges, calib_settings.use_scalers)

        # n_trials = 5
        trial_contribution = Int64[]  # zeros(Int64, n_trials)
        n_trials = 0

        initial_guess = if isfile("./data/ensemble/16071S_initial_guess.dat")
            deserialize("./data/ensemble/16071S_initial_guess.dat")
        else
            # Will error. If restarting calibration, best to temporarily comment
            # out the `population=guess` argument.
            # Sorry, don't have time to make this more flexible.
            nothing
        end

        # Run multiple optimizations with different starting conditions
        # for trial in 1:n_trials
        while sum(trial_contribution) < opt_config.ensemble_members
            n_trials += 1
            @info "Running attempt $(n_trials)"
            @info "Progress so far: $(trial_contribution) [Total: $(sum(trial_contribution))]"
            trial_counter[] = 0

            guess = if !isnothing(initial_guess)
                guess_subset = rand(
                    1:size(initial_guess, 2),
                    floor(Int64, opt_config.population_size * 0.5)
                )
                initial_guess[:, guess_subset]
            else
                nothing
            end

            # Run optimization
            res = bboptimize(
                objective;
                SearchRange=collect(search_ranges),
                MaxSteps=opt_config.max_steps,
                MaxTime=40 * 60,
                Population=guess,
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
    CoralFlow.set_population!(reef_state, optim_best)

    rng = Random.seed!(opt_config.random_seed)
    if calib_settings.use_scalers && length(optim_best) > 16
        n_grps = CoralFlow.n_groups(reef_state)

        # Extract and apply scalers
        scaler_start = 17
        scaler_end = 17 + n_grps - 1
        loc_scalers = optim_best[scaler_start:scaler_end]
        CoralFlow.assign_scalers!(reef_state, loc_scalers)

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

        CoralFlow.run_example!(
            reef_state, env_conditions;
            recruits=Float32(recruitment_proportion),
            self_seed=Float32(self_seeding_proportion),
            rng=rng
        )
    else
        CoralFlow.run_example!(reef_state, env_conditions; rng=rng)
    end

    # Calculate performance metrics for entire time series
    cover = (CoralFlow.coral_cover(reef_state) ./ reef_config.area) * 100.0
    sim = cover[c_sim_indices]
    obs = reef_obs.MEAN_LIVE_CORAL[c_ref_indices]

    calib_metrics = calculate_performance_metrics(sim, obs)

    @info "Performance Metrics:"
    @info "  RMSE: $(round(calib_metrics.rmse; digits=2))%"
    @info "  Pearson: $(round(calib_metrics.pearson; digits=3))"
    @info "  Kendall: $(round(calib_metrics.kendall; digits=3))"
    @info "  Bias (β): $(round(calib_metrics.bias; digits=3))"
    @info "  Variability (α): $(round(calib_metrics.variability; digits=3))"

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
    f_ts = CoralFlow.viz.timeseries(reef_state, env_conditions)
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
        ensemble_res = CoralFlow.run_ensemble!(reef_state, env_conditions, suitable_params)

        # Ensemble timeseries
        f_ensemble = CoralFlow.viz.ensemble_timeseries(
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

    # Calibration comparison plot (with ensemble if available)
    f_calib = plot_calibration_results(
        reef_state, env_conditions, reef_obs, c_sim_indices, c_ref_indices,
        v_sim_idx, v_ref_idx,
        sim_year_range, cover, reef_config.area, calib_metrics,
        joinpath(
            file_paths.figure_dir, "$(reef_config.reef_id)_calibration_comparison.png"
        );
        ensemble_res=ensemble_res
    )

    @info "Calibration complete!"
    return (
        reef_state=reef_state,
        env_conditions=env_conditions,
        results=res,
        best_params=optim_best,
        metrics=calib_metrics,
        tracked_candidates=tracked_candidates,
        tracked_fitnesses=tracked_fitnesses,
        trial_contribution=trial_contribution,
        ensemble_res=ensemble_res
    )
end

ensemble_data_dir = joinpath(OUTPUT_DIR, "ensemble")
mkpath(ensemble_data_dir)

# Run calibration
calibration_output = run_calibration(
    reef_config, file_paths, opt_config, param_bounds, calib_settings
)

output_path = joinpath(ensemble_data_dir, "$(reef_config.reef_id)_ensemble_output.dat")
serialize(output_path, calibration_output)

ensemble_params = hcat(calibration_output.tracked_candidates...)
ensemble_fitnesses = calibration_output.tracked_fitnesses

# Could save found candidates to use as initial guess for later optimization runs
serialize("./data/ensemble/16071S_initial_guess.dat", ensemble_params)

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

median_fitness = median(ensemble_fitnesses)
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
        color=(:gray, 0.1),  # 30% opacity
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

save("$(FIG_DIR)/$(reef_config.reef_id)_calib_param_interactions.png", f; px_per_unit=DPI)
