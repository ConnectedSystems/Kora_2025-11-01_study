using Distributed
using BlackBoxOptim
using Serialization
using DataStructures

ENV["JULIA_DEPOT_PATH"] = "./data/jl_depot"

# Add workers with the current project environment
n_workers = 20  # or use: addprocs(Sys.CPU_THREADS ÷ 2)

if length(workers()) >= n_workers
    rmprocs(workers())
end

# CRITICAL: Activate the same project environment on all workers
addprocs(n_workers; exeflags="--project=@.")

# Load all necessary packages and code on all workers
@everywhere begin
    using BlackBoxOptim
    using CSV
    using DataFrames
    using DataStructures
    using Dates
    using Random
    using Serialization
    
    # Include your modules
    include("common.jl")
    include("_16071S_config.jl")
    
    # Any other includes needed by your objective function
    using CoralFlow
end

"""
    run_single_trial(trial_id, objective, search_ranges, opt_config, 
                     fitness_threshold, ensemble_members, initial_guess)

Run a single optimization trial. Returns (trial_id, res, tracked_candidates, 
tracked_fitnesses, trial_counter).
"""

@everywhere function run_single_trial(
    trial_id::Int,
    objective::Function,
    search_ranges::NamedTuple,
    opt_config::OptimizationConfig,
    fitness_threshold::Float64,
    ensemble_members::Int,
    initial_guess::Union{Nothing,Matrix{Float64}}
)
    @info "Worker $(myid()): Starting trial $(trial_id)"
    
    # Create tracking callback for this trial
    callback, tracked_candidates, tracked_fitnesses, trial_counter = create_tracking_callback(
        fitness_threshold, ensemble_members
    )

    # Set common config
    opts = Dict(
        :SearchRange=>collect(search_ranges),
        :MaxSteps=>opt_config.max_steps,
        :MaxTime=>60*60,
        :PopulationSize=>opt_config.population_size,
        :CallbackInterval=>0.0,
        :CallbackFunction=>callback,
        :TraceInterval=>opt_config.trace_interval
    )
    
    # Generate population guess if initial_guess is provided
    if !isnothing(initial_guess)
        # guess_subset = rand(
        #     1:size(initial_guess, 2),
        #     floor(Int64, opt_config.population_size * 0.8)
        # )

        # Constrain to given bounds
        # (optimization process may produce values that are a little bit out of bounds)
        lb, ub = [first(y) for y in search_ranges], [last(y) for y in search_ranges]
        initial_guess = hcat(map(x -> clamp.(x, lb, ub), eachcol(initial_guess))...)
        opts[:Population] = initial_guess  # [:, guess_subset]
    end

    # Run optimization
    res = bboptimize(
        objective;
        opts...
    )
    
    n_found = trial_counter[]
    @info "Worker $(myid()): Trial $(trial_id) complete - found $(n_found) candidates"
    
    return (
        trial_id=trial_id,
        res=res,
        tracked_candidates=tracked_candidates,
        tracked_fitnesses=tracked_fitnesses,
        n_candidates=n_found
    )
end

"""
    merge_best_candidates!(main_candidates, main_fitnesses, new_candidates, new_fitnesses, max_size)

Merge new candidates into main buffers, keeping only the best solutions by fitness.
"""

function merge_best_candidates!(
    main_candidates::CircularBuffer{Vector{Float64}},
    main_fitnesses::CircularBuffer{Float64},
    new_candidates,
    new_fitnesses,
    max_size::Int
)
    # Combine existing and new candidates
    all_candidates = vcat(collect(main_candidates), collect(new_candidates))
    all_fitnesses = vcat(collect(main_fitnesses), collect(new_fitnesses))
    
    # Sort by fitness (lower is better)
    sorted_idx = sortperm(all_fitnesses)
    
    # Keep only the best max_size candidates
    n_keep = min(length(sorted_idx), max_size)
    best_idx = sorted_idx[1:n_keep]
    
    # Clear and refill buffers with best candidates
    empty!(main_candidates)
    empty!(main_fitnesses)
    
    for idx in best_idx
        push!(main_candidates, all_candidates[idx])
        push!(main_fitnesses, all_fitnesses[idx])
    end
    
    return length(best_idx)
end

"""
    run_calibration_parallel(reef_config, file_paths, opt_config, search_ranges, 
                             calib_settings; max_trials=40)

Run complete ensemble workflow with parallel trials. Runs trials across all available
workers until target number of ensemble members is reached.
"""

@everywhere function run_calibration_parallel(
    reef_config::ReefConfig,
    file_paths::CalibrationDataPaths,
    opt_config::OptimizationConfig,
    search_ranges::NamedTuple,
    calib_settings::CalibrationSettings;
    max_trials::Int=40  # Safety limit to prevent infinite loops
)
    @info "Starting parallel ensemble search for $(reef_config.reef_name) ($(reef_config.reef_id))"
    
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
    ltmp_start_year = Year(Date(reef_obs.SAMPLE_DATE[3])).value
    ltmp_end_year = Year(Date(reef_obs.SAMPLE_DATE[end])).value
    
    benthic_estimate = CSV.read("$(OUTPUT_DIR)/ecorrap_benthic/moore_estimate.csv", DataFrame)
    
    sim_start_year = ltmp_start_year
    sim_end_year = benthic_estimate.year[end]
    
    historic_dhw = CSV.read("$(OUTPUT_DIR)/dhw/DHW_Moore_Reef.csv", DataFrame)
    historic_dhw = historic_dhw[historic_dhw.year .∈ Ref(sim_start_year:sim_end_year), "dhw"]
    
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
        sim_start_year, sim_end_year, calib_settings.start_month
    )
    c_sim_indices, c_ref_indices, c_matched_dates = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
    )
    
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
    
    # Get benthic data used for calibration
    year_span = year.(c_matched_dates)
    sim_benthic_years = [year_span .∈ Ref(benthic_estimate.year)][1]
    aligned_years = year_span[sim_benthic_years]
    benthic_data = benthic_estimate[benthic_estimate.year .∈ aligned_years, :]
    
    if all(isfile.(result_files))
        @info "Loading existing ensemble results..."
        res, optim_best, tracked_candidates, tracked_fitnesses, trial_contribution = load_calibration_results(
            ensemble_dir, reef_config.reef_id
        )
        opt_results = []  # Not stored when loading
        total_trials = length(trial_contribution)
    else
        @info "Running parallel optimization..."
        
        # Create objective function
        objective = create_objective_function(
            reef_state, env_conditions, reef_obs.MEAN_LIVE_CORAL,
            c_sim_indices, c_ref_indices, reef_config.area, sim_benthic_years, benthic_data,
            opt_config.random_seed, calib_settings.use_scalers
        )
        
        # Load or create initial guess
        initial_guess = if isfile("./data/ensemble/16071S_initial_guess.dat")
            deserialize("./data/ensemble/16071S_initial_guess.dat")
        else
            nothing
        end
        
        # Initialize tracking with CircularBuffers from the start
        tracked_candidates = CircularBuffer{Vector{Float64}}(opt_config.ensemble_members)
        tracked_fitnesses = CircularBuffer{Float64}(opt_config.ensemble_members)
        trial_contribution = Int64[]
        opt_results = []
        
        # Run trials until we have enough ensemble members
        trial_id = 1
        n_workers_available = length(workers())
        
        while length(tracked_candidates) < opt_config.ensemble_members && trial_id <= max_trials
            # Determine how many trials to run in this iteration
            n_trials_to_run = min(n_workers_available, max_trials - trial_id + 1)
            trial_ids = trial_id:(trial_id + n_trials_to_run - 1)
            
            @info "=" ^ 60
            @info "Running trials $(first(trial_ids)):$(last(trial_ids)) in parallel"
            @info "Current progress: $(length(tracked_candidates)) / $(opt_config.ensemble_members) members"
            
            # Run trials in parallel using pmap
            trial_results = pmap(id -> run_single_trial(
                id,
                objective,
                search_ranges,
                opt_config,
                opt_config.fitness_threshold,
                opt_config.ensemble_members,
                initial_guess
            ), trial_ids)
            
            # Aggregate results - merge new candidates into main buffers, keeping best
            for result in trial_results
                push!(opt_results, result.res)
                push!(trial_contribution, result.n_candidates)
                
                # Merge candidates from this trial, keeping only the best
                if result.n_candidates > 0
                    merge_best_candidates!(
                        tracked_candidates,
                        tracked_fitnesses,
                        result.tracked_candidates,
                        result.tracked_fitnesses,
                        opt_config.ensemble_members
                    )
                end
            end
            
            trial_id += n_trials_to_run
            
            @info "Iteration complete:"
            @info "  Total trials run: $(trial_id - 1)"
            @info "  Total candidates found: $(length(tracked_candidates))"
            @info "  Last trial contributions: $(trial_contribution[end-n_trials_to_run+1:end])"
            
            # Check completion
            if length(tracked_candidates) >= opt_config.ensemble_members
                @info "Target reached! Found $(length(tracked_candidates)) candidates (need $(opt_config.ensemble_members))"
                break
            end
        end
        
        total_trials = trial_id - 1
        
        if length(tracked_candidates) < opt_config.ensemble_members
            @warn "Reached maximum trials ($(max_trials)) without finding enough candidates"
            @warn "Found $(length(tracked_candidates)) / $(opt_config.ensemble_members) candidates"
        end
        
        # Get best overall candidate
        optim_best = if length(tracked_candidates) > 0
            best_idx = argmin(collect(tracked_fitnesses))
            collect(tracked_candidates)[best_idx]
        else
            # Get best from any trial if no ensemble
            best_candidate(opt_results[1])
        end
        
        # Convert probability values to realized values
        optim_best[2:6] = gamma_to_dirichlet(optim_best[2:6])
        for c in tracked_candidates
            c[2:6] .= gamma_to_dirichlet(c[2:6])
        end
        
        # Save results
        save_calibration_results(
            ensemble_dir, reef_config.reef_id, opt_results[end], optim_best,
            tracked_candidates, tracked_fitnesses, trial_contribution,
            opt_config.fitness_threshold
        )
        
        res = opt_results[end]  # For compatibility with rest of code
    end
    
    if !isempty(opt_results)
        @info "Best fitness across all trials: $(minimum(best_fitness.(opt_results)))"
    end
    
    # Run model with best parameters
    @info "Running model with best parameters..."
    CoralFlow.set_population!(reef_state, optim_best)
    
    rng = Random.seed!(opt_config.random_seed)
    if calib_settings.use_scalers && length(optim_best) > 16
        n_grps = CoralFlow.n_groups(reef_state)
        
        scaler_start = 17
        scaler_end = 17 + n_grps - 1
        loc_scalers = optim_best[scaler_start:scaler_end]
        CoralFlow.assign_scalers!(reef_state, loc_scalers)
        
        recruitment_proportion = optim_best[scaler_end + 1]
        self_seeding_proportion = optim_best[scaler_end + 2]
        
        dhw_tols = extract_dhw_tolerances(optim_best, n_grps)
        if !isnothing(dhw_tols)
            dhw_means, dhw_stds = dhw_tols
            for grp in 1:n_grps
                reef_state.wild_dhw_tolerances[1, :, grp, At(:mean)] .= dhw_means[grp]
                reef_state.wild_dhw_tolerances[:, :, grp, At(:stdev)] .= dhw_stds[grp]
            end
        end
        
        CoralFlow.run_model!(
            reef_state, env_conditions;
            recruits=Float32(recruitment_proportion),
            self_seed=Float32(self_seeding_proportion),
            rng=rng
        )
    else
        CoralFlow.run_model!(reef_state, env_conditions; rng=rng)
    end
    
    # Calculate metrics and generate visualizations
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
    
    v_sim_indices, v_ref_indices, _ = find_closest_dates(
        collect(Date.(sim_year_range)),
        Date.(reef_obs.SAMPLE_DATE);
        max_days=calib_settings.date_match_max_days
    )
    
    v_sim_idx = [s for s in v_sim_indices if s ∉ c_sim_indices]
    v_ref_idx = [s for s in v_ref_indices if s ∉ c_ref_indices]
    
    @info "Generating visualizations..."
    mkpath(file_paths.figure_dir)
    
    f_ts = CoralFlow.viz.timeseries(reef_state, env_conditions)
    save(
        joinpath(file_paths.figure_dir, "$(reef_config.reef_id)_calibrated_timeseries.png"),
        f_ts;
        px_per_unit=DPI
    )
    
    ensemble_res = nothing
    if !isempty(tracked_candidates)
        @info "Running ensemble with $(length(tracked_candidates)) candidates..."
        suitable_params = hcat(collect(tracked_candidates)...)
        ensemble_res = CoralFlow.run_ensemble!(reef_state, env_conditions, suitable_params)
        
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
    
    f_calib = plot_calibration_results(
        reef_state, env_conditions, reef_obs, c_sim_indices, c_ref_indices,
        sim_year_range, cover, reef_config.area, benthic_estimate,
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
        tracked_candidates=collect(tracked_candidates),  # Convert to Vector for return
        tracked_fitnesses=collect(tracked_fitnesses),
        trial_contribution=trial_contribution,
        ensemble_res=ensemble_res,
        n_trials=total_trials
    )
end

# Example usage
ensemble_data_dir = joinpath(OUTPUT_DIR, "ensemble")
mkpath(ensemble_data_dir)

# Run parallel calibration
calibration_output = run_calibration_parallel(
    reef_config, file_paths, opt_config, param_bounds, calib_settings;
    max_trials=5_000_000  # Adjust based on computational budget
)

# Save outputs
output_path = joinpath(ensemble_data_dir, "$(reef_config.reef_id)_ensemble_output.dat")
serialize(output_path, calibration_output)

ensemble_params = hcat(calibration_output.tracked_candidates...)
ensemble_fitnesses = calibration_output.tracked_fitnesses

serialize("./data/ensemble/16071S_initial_guess.dat", ensemble_params)

# Parameter identifiability analysis
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

# Visualization
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

save("$(FIG_DIR)/$(reef_config.reef_id)_calib_param_interactions.png", f; px_per_unit=DPI)
