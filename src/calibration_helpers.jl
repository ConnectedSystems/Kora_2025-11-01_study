import Kora: ReefState
import DataStructures: CircularBuffer, capacity

using LaTeXStrings

"""
    load_reef_observations(parquet_file, reef_id)

Load reef observations from MantaTow parquet file.
"""
function load_reef_observations(parquet_file::String, reef_id::String)
    ds = Parquet2.Dataset(parquet_file)
    df = DataFrame(ds; copycols=false)
    gdf = groupby(df, :REEF_ID)
    close(ds)

    reef_obs = gdf[(reef_id,)][:, [:SAMPLE_DATE, :MEAN_LIVE_CORAL]]
    reef_obs.MEAN_LIVE_CORAL .= Float64.(reef_obs.MEAN_LIVE_CORAL)

    return reef_obs
end

"""
    get_reef_uid(canonical_reefs_file, reef_id)

Get the unique reef ID from canonical reefs file.
"""
function get_reef_uid(canonical_reefs_file::String, reef_id::String)
    canonical_reefs = GDF.read(canonical_reefs_file)
    ltmp_reef_ids = replace.(canonical_reefs.LTMP_ID, "-" => "")
    reef_uid = canonical_reefs[occursin.(reef_id, ltmp_reef_ids), :UNIQUE_ID]
    return first(reef_uid)
end

"""
    load_historical_dhw(dhw_file, reef_uid, start_year, end_year; scenario=1)

Load historical DHW data for a single location over the specified time period.

Reads the `dhw_scens` variable from a DHW scenario NetCDF, selecting the location
matching `reef_uid` and all timesteps in `start_year:end_year` (inclusive).

# Returns
`Vector` of DHW values, one per timestep.
"""
function load_historical_dhw(
    dhw_file::String,
    reef_uid,
    start_year::Int,
    end_year::Int;
    scenario::Int=1
)
    return NCDataset(dhw_file, "r") do ds
        dhw_scens = ds["dhw_scens"]

        locs = ds["locs"][:]
        loc_idx = findfirst(==(reef_uid), locs)
        if isnothing(loc_idx)
            throw(ArgumentError("Location $(reef_uid) not found in $(dhw_file)"))
        end

        timesteps = ds["timesteps"][:]
        year_idx = findall(start_year .<= timesteps .<= end_year)
        if isempty(year_idx)
            throw(
                ArgumentError(
                    "No timesteps within $(start_year):$(end_year) found in $(dhw_file)"
                )
            )
        end

        # Index by dimension name so layout changes in the source file are tolerated
        selection = Dict(
            "timesteps" => year_idx, "locs" => loc_idx, "scenarios" => scenario
        )
        return vec(dhw_scens[getindex.(Ref(selection), NCDatasets.dimnames(dhw_scens))...])
    end
end

"""
    create_simulation_dates(start_year, end_year, start_month=3)

Create date range for simulation starting in specified month.
"""
function create_simulation_dates(start_year::Int, end_year::Int, start_month::Int=3)
    start_date = Date("$(start_year)-$(lpad(start_month, 2, '0'))-01")
    end_date = Date("$(end_year)-$(lpad(start_month, 2, '0'))-01")
    return start_date:Year(1):end_date
end

"""
    find_closest_dates(target_dates, reference_dates; max_days=120)

Find closest reference date for each target date. Ignored if a date within
max_days is not found.

# Returns
- Indices for simulation for comparison
- Indices of reference dates to use for comparison
- Reference dates used for comparison
"""
function find_closest_dates(
    target_dates::AbstractVector{Date},
    reference_dates::AbstractVector{Date};
    max_days::Int=120
)
    closest_ref_indices = Vector{Union{Int,Missing}}(missing, length(target_dates))
    closest_ref_dates = Vector{Union{Date,Missing}}(missing, length(target_dates))

    for (i, target) in enumerate(target_dates)
        diffs = abs.(reference_dates .- target)
        min_idx = argmin(diffs)

        # Find observations within n days of target date
        # otherwise, leave it as `missing`
        if diffs[min_idx] > Day(max_days)
            continue
        end

        closest_ref_indices[i] = min_idx
        closest_ref_dates[i] = reference_dates[min_idx]
    end

    return (
        findall(.!ismissing.(closest_ref_indices)),
        collect.([
            skipmissing(closest_ref_indices),
            skipmissing(closest_ref_dates)
        ])...
    )
end

"""
    exclude_years_from_indices(sim_indices::Vector{Int}, ref_indices::Vector{Int}, matched_dates::Vector{Date}, years_to_exclude::Vector{Int})

Remove specified years from calibration indices.

# Arguments
- `sim_indices`: Indices of simulation timepoints to compare
- `ref_indices`: Indices of observation timepoints to compare
- `matched_dates`: Matched observation dates
- `years_to_exclude`: Years to exclude from calibration
"""
function exclude_years_from_indices(
    sim_indices::Vector{Int},
    ref_indices::Vector{Int},
    matched_dates::Vector{Date},
    years_to_exclude::Vector{Int}
)
    if isempty(years_to_exclude)
        return sim_indices, ref_indices, matched_dates
    end

    exclude_years = [Year(m).value in years_to_exclude for m in matched_dates]

    return (
        sim_indices[Not(exclude_years)],
        ref_indices[Not(exclude_years)],
        matched_dates[Not(exclude_years)]
    )
end

"""
    calculate_slope_penalty(covers; window_size=10, threshold=0.2)

Calculate penalty for implausibly smooth trajectories.
"""
function calculate_slope_penalty(
    covers::Matrix{Float32};
    window_size::Int=10,
    threshold::Float32=0.2f0
)
    tf = size(covers, 1)
    slope_penalty = 0.0

    for g_id in axes(covers, 2)
        ts = covers[:, g_id]
        group_max_penalty = 0.0

        for i in 1:(tf - window_size + 1)
            w = ts[i:(i + window_size - 1)]
            slopes = diff(w)
            mean_slope = mean(abs.(slopes))
            slope_std = std(slopes)

            cv = mean_slope > 1e-6 ? slope_std / mean_slope : 0.0f0
            penalty = cv < threshold ? (threshold - cv) : 0.0f0
            group_max_penalty = max(group_max_penalty, penalty)
        end

        slope_penalty += group_max_penalty
    end

    return slope_penalty
end

"""
    gamma_to_dirichlet(p; G=Gamma(1))

Gamma-Dirichlet transformation.

Use Gamma distribution trick to convert values `p` into those that sum to 1.

Treat vector `p` as probability levels (0 - 1), used to compute corresponding
quantiles from a Gamma distribution. These are normalized so that the transformation sums
to 1.

Taking the quantiles from a Gamma distribution with α:=1 creates a uniform Dirichlet sample.

# Arguments
- `p` : Vector of probability levels
- `G` : Gamma distribution (default: Gamma(α=1))

# Returns
Vector of values that sum to 1.
"""
function gamma_to_dirichlet(p; G=Gamma(1))
    gamma_samples = quantile.(G, p)
    return gamma_samples ./ sum(gamma_samples)
end

"""
    dirichlet_to_latent(p)

Inverse of `gamma_to_dirichlet` (with the default `Gamma(1)`) for proportions `p` that
sum to 1: `gamma_to_dirichlet(dirichlet_to_latent(p)) ≈ p`. The latent representation is
not unique (any positive rescaling of the Gamma quantiles gives the same proportions);
this returns the one whose quantiles equal `p` itself.
"""
dirichlet_to_latent(p) = 1 .- exp.(-p)

"""
    candidate_rng(x, seed)

Random number generator for evaluating parameter vector `x`, derived from `x` and `seed`.
The same candidate and seed always yield the same random stream, independent of thread,
evaluation order, or which candidates were evaluated before it, while distinct candidates
receive independent streams.
"""
candidate_rng(x, seed::Integer) = Xoshiro(hash(x, hash(seed)))

"""
    create_objective_function(reef_state, env_conditions, historic_obs, sim_indices,
                              ref_indices, area, seed, use_scalers;
                              transform_proportions=true)

Create objective function for calibration optimization.

With `transform_proportions=true` (calibration), `x[2:6]` are latent values converted
to proportions via `gamma_to_dirichlet`. Pass `false` when `x[2:6]` already hold
realized proportions (e.g., stored ensemble candidates or pre-transformed SA samples),
so they are used as-is rather than transformed a second time.

Randomness comes from `candidate_rng(x, seed)`, so the objective is a deterministic
function of `x` for a given `seed`.
"""
function create_objective_function(
    reef_state::ReefState,
    env_conditions::AbstractDimArray,
    historic_obs::Vector{Float64},
    sim_indices::Vector{Int},
    ref_indices::Vector{Int},
    area::Float32,
    benthic_aligned_years::Union{BitVector,Nothing},
    benthic::Union{DataFrame,Nothing},
    seed::Int,
    use_scalers::Bool;
    transform_proportions::Bool=true
)
    obs = historic_obs[ref_indices] / 100.0
    obs_μ, obs_σ = mean_and_std(obs)

    function objective(x)
        # Use Gamma distribution trick to convert candidate values to those that sum to 1.
        # A copy is necessary as values will be adjusted
        x = if transform_proportions
            vcat(x[1], gamma_to_dirichlet(x[2:6]), x[7:end])
        else
            copy(x)
        end

        # Alternate approach: penalize invalid combinations
        # (computationally wasteful, the gamma transform is better)
        # if !(sum(x[2:6]) ≈ 1.0)
        #     return 1e10 * abs(1.0 - sum(x[2:6]))
        # end

        rng = candidate_rng(x, seed)
        Kora.set_population!(reef_state, x; rng=rng)

        if use_scalers
            n_grps = Kora.n_groups(reef_state)
            scaler_end = 17 + n_grps - 1
            loc_scalers = x[17:scaler_end]
            Kora.assign_scalers!(reef_state, loc_scalers)

            recruitment_proportion = x[end - 1]
            self_seeding_proportion = x[end]

            Kora.run_model!(
                reef_state,
                env_conditions;
                recruits=Float32(recruitment_proportion),
                self_seed=Float32(self_seeding_proportion),
                rng=rng
            )
        else
            Kora.run_model!(reef_state, env_conditions; rng=rng)
        end

        cover = Kora.coral_cover(reef_state)
        cover ./= area
        sims = cover[sim_indices]

        # Return a large penalty for degenerate trajectories (NaN/Inf cover or
        # zero-variance sims, which would produce NaN from pearson/α/β).
        if any(!isfinite, sims) || std(sims) == 0.0
            return 1e6
        end

        # Error metrics
        init_mae = abs(sims[1] - obs[1])
        end_mae = abs(sims[end] - obs[end])
        We = (init_mae + end_mae)
        # rmse = Kora.RMSE(sims, obs)
        pearson = (1.0 - Kora.pearson(sims, obs)) / 2.0

        sim_μ, std_hat = mean_and_std(sims)

        # Bias metric
        u = abs((sim_μ / obs_μ) - 1.0)
        β = u / (1 + u)

        # Variability ratio
        _v = abs((std_hat / obs_σ) - 1.0)
        α = (_v / (1 + _v))

        # Trajectory plausibility penalty
        # slope_penalty = calculate_slope_penalty(covers)

        # Penalize trajectories where a group falls below 1% more than 3 times
        covers = Kora.group_cover_timeseries(reef_state)
        low_cover_penalty = sum(count(covers .< 1.0; dims=1) .> 3)

        # Alignment with benthic observations
        rank_score = 0.0
        if !isnothing(benthic_aligned_years)
            benthic_sim = @view covers[sim_indices, :][benthic_aligned_years, :]

            if any(benthic_sim .== 0.0)
                rank_score = 2.0 * (count(benthic_sim .== 0.0))
            else
                # Pre-allocate vectors for correlation calculation
                s_buffer = Vector{Float64}(undef, size(benthic_sim, 2))
                o_buffer = Vector{Float64}(undef, size(benthic, 2) - 1)

                # Assess for each year where benthic obs are available
                for (i, r) in enumerate(eachrow(benthic))
                    s_buffer .= @view benthic_sim[i, :]
                    o_buffer .= collect(r[2:end])
                    log_r = std(s_buffer) == 0.0 ? 2.0 : -cor(log.(s_buffer), log.(o_buffer)) + 1.0
                    rank_score += log_r
                end

                rank_score /= nrow(benthic)  # average of years
                # rank_score *= 0.25  # place less weight on rank order
            end
        end

        return (α + β + pearson) + We + low_cover_penalty + rank_score
    end

    return objective
end

"""
    create_objective_pool(reef_state, env_conditions, historic_obs, sim_indices,
                          ref_indices, area, benthic_aligned_years, benthic,
                          seed, use_scalers; transform_proportions=true)

Create a thread-safe pool of objective functions for parallel sensitivity analysis.
Backed by `Threads.OncePerThread`, so each thread lazily builds — on first use — its
own `create_objective_function` closure over an independent `copy` of `reef_state`.
All threads share `seed`, and randomness comes from `candidate_rng(x, seed)`, so results
do not depend on thread assignment or evaluation order. Safe for use with `Threads.@threads`.

# Arguments
All arguments are forwarded to `create_objective_function`.

# Returns
A callable `Threads.OncePerThread`; call it with no arguments to get the current
thread's objective function, e.g. `pool()(x)`.
"""
function create_objective_pool(
    reef_state::ReefState,
    env_conditions::AbstractDimArray,
    historic_obs::Vector{Float64},
    sim_indices::Vector{Int},
    ref_indices::Vector{Int},
    area::Float32,
    benthic_aligned_years::Union{BitVector,Nothing},
    benthic::Union{DataFrame,Nothing},
    seed::Int,
    use_scalers::Bool;
    transform_proportions::Bool=true
)
    return Threads.OncePerThread{Function}() do
        create_objective_function(
            copy(reef_state), env_conditions, historic_obs,
            sim_indices, ref_indices, area,
            benthic_aligned_years, benthic,
            seed, use_scalers;
            transform_proportions=transform_proportions
        )
    end
end

"""
    create_timeseries_function(reef_state, env_conditions, area, seed, use_scalers;
                               transform_proportions=true)

Create a function that runs the model and returns the full simulated coral cover
trajectory rather than aggregating to a scalar fitness score. Shares the same
parameter layout and `transform_proportions` handling as `create_objective_function`.
Stored `ensemble_params` columns hold realized proportions, so pass
`transform_proportions=false` when running them.

Thread-safe: `reef_state` is mutated in-place by `run_model!`, so each thread
lazily gets its own `copy` via `Threads.OncePerThread`. By default randomness comes from
`candidate_rng(x, seed)`, so a candidate gives the same trajectory as its evaluation by
`create_objective_function` with the same `seed`, and total and by-group runs match.
Pass `rng` to use a given random stream instead (e.g. to pair runs of one sample under two
scenarios; see `recruitment_paired_runs`).

# Returns
A function `(x; return_by_group=false, rng=nothing) → cover`, where `cover` is total cover
(`Vector{Float32}`, length `n_timesteps`) or per-group cover (`Matrix{Float32}`,
`n_timesteps × n_groups`) if `return_by_group=true`.
"""
function create_timeseries_function(
    reef_state::ReefState,
    env_conditions::AbstractDimArray,
    area::Float32,
    seed::Int,
    use_scalers::Bool;
    transform_proportions::Bool=true
)
    state_pool = Threads.OncePerThread{typeof(reef_state)}(() -> copy(reef_state))

    function timeseries(
        x; return_by_group::Bool=false, rng::Union{Nothing,AbstractRNG}=nothing
    )
        rs = state_pool()

        x = if transform_proportions
            vcat(x[1], gamma_to_dirichlet(x[2:6]), x[7:end])
        else
            copy(x)
        end
        rng = isnothing(rng) ? candidate_rng(x, seed) : rng

        Kora.set_population!(rs, x; rng=rng)

        if use_scalers
            n_grps = Kora.n_groups(rs)
            scaler_end = 17 + n_grps - 1
            loc_scalers = x[17:scaler_end]
            Kora.assign_scalers!(rs, loc_scalers)
            recruitment_proportion = x[end - 1]
            self_seeding_proportion = x[end]
            Kora.run_model!(
                rs, env_conditions;
                recruits=Float32(recruitment_proportion),
                self_seed=Float32(self_seeding_proportion),
                rng=rng
            )
        else
            Kora.run_model!(rs, env_conditions; rng=rng)
        end

        if return_by_group
            # (n_timesteps × n_groups) — enables per-group temporal PAWN
            return Kora.group_cover_timeseries(rs) ./ area
        else
            cover = Kora.coral_cover(rs)
            return cover ./ area
        end
    end

    return timeseries
end

"""
    recruitment_paired_runs(timeseries, X, seed) → (with, without)

Run every sample (row of `X`) twice through `timeseries` (from `create_timeseries_function`
with `use_scalers=true`): once as given, and once with both recruitment parameters (external
recruitment `x[end - 1]` and self-seeding `x[end]`) set to zero.

Both runs of sample `i` use `Xoshiro(hash(i, hash(seed)))`, so they share the initial
population and random stream; their difference is the cover attributable to recruitment
rather than to model noise. (The streams can still diverge once recruitment changes the
sequence of random draws.) Returns two `n_samples × n_timesteps` matrices of total cover.
"""
function recruitment_paired_runs(timeseries, X::AbstractMatrix, seed::Integer)
    n_samples = size(X, 1)
    n_steps = length(timeseries(X[1, :]))
    with_recruitment = Matrix{Float32}(undef, n_samples, n_steps)
    without_recruitment = similar(with_recruitment)

    Threads.@threads :dynamic for i in 1:n_samples
        x = X[i, :]
        with_recruitment[i, :] = timeseries(x; rng=Xoshiro(hash(i, hash(seed))))

        x[end - 1] = 0.0
        x[end] = 0.0
        without_recruitment[i, :] = timeseries(x; rng=Xoshiro(hash(i, hash(seed))))
    end

    return with_recruitment, without_recruitment
end

"""
    create_tracking_callback(fitness_threshold)

Create callback function to track good optimization candidates.
"""
function create_tracking_callback(fitness_threshold::Float64, n_members::Int64)
    tracked_candidates = CircularBuffer{Vector{Float64}}(n_members)
    tracked_fitnesses = CircularBuffer{Float64}(n_members)
    trial_counter = Ref(0)
    candidate_progress = CircularBuffer{Int64}(10000)

    function callback(oc)
        # Track num_better (which will be 0 if no change)
        push!(candidate_progress, oc.num_better)

        elapsed_time = (time() - oc.start_time)

        # No progress early on
        if (30 < elapsed_time < 200) && (length(candidate_progress) >= 1000)
            no_progress = all(candidate_progress .== candidate_progress[end])
            if no_progress
                BlackBoxOptim.shutdown!(oc)
                println("Early stopping: No progress in last 1000 iterations")
                return true
            end
        end

        archive = try
            oc.optimizer.population
        catch
            # For separable_nes method
            (;
                fitness=fitness.(oc.optimizer.candidates),
                individuals=hcat(getfield.(oc.optimizer.candidates, :params)...)
            )
        end

        for (i, fitness) in enumerate(archive.fitness)
            if !(fitness < fitness_threshold)
                continue
            end

            candidate = archive.individuals[:, i]

            if isempty(tracked_candidates)
                # Add and move on
                push!(tracked_candidates, copy(candidate))
                push!(tracked_fitnesses, fitness)

                println("Added first candidate! Fitness: $(fitness)")
                trial_counter[] += 1

                continue
            end

            # Check for duplicates
            is_duplicate = any(tracked_candidates) do existing
                all(abs.(existing .- candidate) .< 1e-8)
            end

            if !is_duplicate
                n_tracked = length(tracked_candidates)
                if n_tracked == capacity(tracked_candidates)
                    # Exit once ensemble collated
                    # BlackBoxOptim.shutdown!(oc)

                    # Has to beat at least the bottom 20%
                    if fitness < quantile(tracked_fitnesses, 0.8)
                        # Replace the worst scoring candidate
                        _, idx = findmax(tracked_fitnesses)
                        tracked_candidates[idx] = copy(candidate)
                        tracked_fitnesses[idx] = fitness
                        trial_counter[] += 1
                    else
                        continue
                    end
                else
                    push!(tracked_candidates, copy(candidate))
                    push!(tracked_fitnesses, fitness)
                    n_tracked += 1
                    trial_counter[] += 1
                end

                if (n_tracked < 10) || trial_counter[] % 10 == 0
                    println(
                        "Tracked $(n_tracked) candidates (added: $(trial_counter[])) " *
                        "(fitness < $fitness_threshold, mean fitness: $(mean(tracked_fitnesses)))"
                    )
                end
            end
        end

        # After 20 minutes check fitness threshold or candidate count
        n_tracked = trial_counter[]
        if elapsed_time >= 600  # 10 minutes
            # Check for progress and quit if stuck in local optima
            no_progress = all(candidate_progress .== candidate_progress[end])
            if no_progress
                BlackBoxOptim.shutdown!(oc)
                return true
            end

            # Otherwise check average fitness
            avg_fitness = mean(oc.optimizer.population.fitness)
            if avg_fitness < 0.15
                println(
                    "Stopping: Average fitness ($(round(avg_fitness, digits=4))) below threshold"
                )
                BlackBoxOptim.shutdown!(oc)
                return true
            end
        end

        if n_tracked >= 20
            println("Stopping: Found $(n_tracked) candidate solutions")
            BlackBoxOptim.shutdown!(oc)
            return true
        end

        return false
    end

    return callback, tracked_candidates, tracked_fitnesses, trial_counter
end

"""
    create_stagnation_callback()

Create a callback that performs early stopping on optimizer stagnation only (no
candidate tracking). Mirrors the stagnation-detection logic in
[`create_tracking_callback`](@ref); used together with
[`harvest_final_population`](@ref), which takes a single end-of-trial snapshot of
the best candidates instead of accumulating them across every iteration.
"""
function create_stagnation_callback()
    candidate_progress = CircularBuffer{Int64}(10000)

    function callback(oc)
        push!(candidate_progress, oc.num_better)

        elapsed_time = (time() - oc.start_time)

        if (30 < elapsed_time < 200) && (length(candidate_progress) >= 1000)
            no_progress = all(candidate_progress .== candidate_progress[end])
            if no_progress
                BlackBoxOptim.shutdown!(oc)
                println("Early stopping: No progress in last 1000 iterations")
                return true
            end
        end

        if elapsed_time >= 600  # 10 minutes
            no_progress = all(candidate_progress .== candidate_progress[end])
            if no_progress
                BlackBoxOptim.shutdown!(oc)
                return true
            end

            avg_fitness = mean(oc.optimizer.population.fitness)
            if avg_fitness < 0.15
                println(
                    "Stopping: Average fitness ($(round(avg_fitness, digits=4))) below threshold"
                )
                BlackBoxOptim.shutdown!(oc)
                return true
            end
        end

        return false
    end

    return callback
end

"""
    harvest_final_population(res, fitness_threshold; n_best=20)

Take a single end-of-trial snapshot of `res`'s final population (as opposed to
[`create_tracking_callback`](@ref)'s per-iteration accumulation across the whole
run), returning the best `n_best` distinct candidates below `fitness_threshold`.

Only distinct parameter sets are returned. A converged population can hold many identical
copies of one individual, and with a deterministic objective each copy has the same
fitness; keeping them would add duplicate members to the ensemble.

Falls back to a `(fitness, individuals)` reconstruction for optimizer methods whose
`population(res)` does not expose `.fitness`/`.individuals` directly (mirrors the
equivalent fallback in [`create_tracking_callback`](@ref)).
"""
function harvest_final_population(res, fitness_threshold::Float64; n_best::Int=20)
    pop = population(res)
    fits, indivs = try
        pop.fitness, pop.individuals
    catch
        fitness.(pop), hcat(getfield.(pop, :params)...)
    end

    below = findall(f -> f < fitness_threshold, fits)
    below = below[sortperm(fits[below])]
    below = below[unique(i -> indivs[:, below[i]], eachindex(below))]
    keep = below[1:min(n_best, length(below))]

    return [indivs[:, i] for i in keep], fits[keep]
end

"""
    build_search_ranges(ranges, use_scalers)

Build search range vector for optimization.

Accepts either a `SearchRanges` struct or a `NamedTuple` with the same fields.
"""
function build_search_ranges(ranges::SearchRanges, use_scalers::Bool)
    base_ranges = vcat(
        [ranges.density],
        fill(ranges.group_proportion, 5),
        fill(ranges.size_mean, 5),
        fill(ranges.size_std, 5)
    )

    if use_scalers
        return vcat(
            base_ranges,
            fill(ranges.scalers, 5),
            [ranges.recruitment],
            [ranges.self_seeding]
        )
    end

    return base_ranges
end

"""
    extract_dhw_tolerances(optim_best, n_groups)

Extract calibrated DHW tolerance means and standard deviations from optimization results.

# Returns
- `Tuple{Vector{Float32}, Vector{Float32}}`: (dhw_means, dhw_stds) for each group
"""
function extract_dhw_tolerances(optim_best::Vector{Float64}, n_groups::Int)
    # Check if we have DHW tolerance parameters
    # Base parameters: 16, Scalers: 5, Recruitment/Self-seed: 2, DHW: 10
    expected_length = 16 + n_groups + 2 + (2 * n_groups)

    if length(optim_best) < expected_length
        # No DHW tolerances in parameters
        return nothing
    end

    # DHW means start after: base(16) + scalers(5) + recruitment/self-seed(2)
    dhw_mean_start = 16 + n_groups + 3
    dhw_mean_end = dhw_mean_start + n_groups - 1
    dhw_means = Float32.(optim_best[dhw_mean_start:dhw_mean_end])

    dhw_std_start = dhw_mean_end + 1
    dhw_std_end = dhw_std_start + n_groups - 1
    dhw_stds = Float32.(optim_best[dhw_std_start:dhw_std_end])

    return (dhw_means, dhw_stds)
end

# """
#     extract_calibration_parameters(x::Vector{Float64}, n_groups::Int)

# Extract all calibration parameters from optimization parameter vector.

# # Arguments
# - `x::Vector{Float64}`: Parameter vector from optimization
# - `n_groups::Int`: Number of functional groups (typically 5)

# # Returns
# Named tuple with:
# - `density::Float64`: Population density
# - `group_proportions::Vector{Float64}`: Group proportions (5 values)
# - `size_means::Vector{Float64}`: Size distribution means (5 values)
# - `size_stds::Vector{Float64}`: Size distribution stdevs (5 values)
# - `scalers::Union{Vector{Float64}, Nothing}`: Growth/survival scalers (5 values or nothing)
# - `recruitment::Union{Float64, Nothing}`: Recruitment proportion (or nothing)
# - `self_seeding::Union{Float64, Nothing}`: Self-seeding proportion (or nothing)
# - `dhw_means::Union{Vector{Float64}, Nothing}`: DHW tolerance means (5 values or nothing)
# - `dhw_stds::Union{Vector{Float64}, Nothing}`: DHW tolerance stdevs (5 values or nothing)
# """
# function extract_calibration_parameters(x::Vector{Float64}, n_groups::Int)
#     # Base parameters (always present)
#     density = x[1]
#     group_proportions = x[2:(1 + n_groups)]
#     size_means = x[(2 + n_groups):(1 + 2 * n_groups)]
#     size_stds = x[(2 + 2 * n_groups):(1 + 3 * n_groups)]

#     base_length = 1 + 3 * n_groups  # 16 for n_groups=5

#     # Check if extended parameters are present
#     if length(x) <= base_length
#         return (
#             density=density,
#             group_proportions=group_proportions,
#             size_means=size_means,
#             size_stds=size_stds,
#             scalers=nothing,
#             recruitment=nothing,
#             self_seeding=nothing,
#             dhw_means=nothing,
#             dhw_stds=nothing
#         )
#     end

#     # Extract scalers and recruitment parameters
#     scaler_start = base_length + 1
#     scaler_end = scaler_start + n_groups - 1
#     scalers = x[scaler_start:scaler_end]

#     recruitment = x[scaler_end + 1]
#     self_seeding = x[scaler_end + 2]

#     # Check if DHW parameters are present
#     dhw_start = scaler_end + 3
#     expected_full_length = scaler_end + 2 + 2 * n_groups

#     if length(x) < expected_full_length
#         return (
#             density=density,
#             group_proportions=group_proportions,
#             size_means=size_means,
#             size_stds=size_stds,
#             scalers=scalers,
#             recruitment=recruitment,
#             self_seeding=self_seeding,
#             dhw_means=nothing,
#             dhw_stds=nothing
#         )
#     end

#     # Extract DHW parameters
#     dhw_mean_end = dhw_start + n_groups - 1
#     dhw_means = x[dhw_start:dhw_mean_end]

#     dhw_std_start = dhw_mean_end + 1
#     dhw_std_end = dhw_std_start + n_groups - 1
#     dhw_stds = x[dhw_std_start:dhw_std_end]

#     return (
#         density=density,
#         group_proportions=group_proportions,
#         size_means=size_means,
#         size_stds=size_stds,
#         scalers=scalers,
#         recruitment=recruitment,
#         self_seeding=self_seeding,
#         dhw_means=dhw_means,
#         dhw_stds=dhw_stds
#     )
# end

"""
    save_calibration_results(output_dir, reef_id, res, optim_best,
                             tracked_candidates, tracked_fitnesses, fitness_threshold)

Save calibration results to disk.
"""
function save_calibration_results(
    output_dir::String,
    reef_id::String,
    res,
    optim_best::Vector{Float64},
    tracked_candidates::CircularBuffer{Vector{Float64}},
    tracked_fitnesses::CircularBuffer{Float64},
    trial_contribution::Vector{Int64},
    fitness_threshold::Float64
)
    mkpath(output_dir)

    prefix = joinpath(output_dir, reef_id)

    # BlackBoxOptim's optimiser state has no array representation, so it stays on
    # Julia serialize; everything else round-trips through HDF5.
    serialize("$(prefix)_optim_state.dat", res)
    save_result("$(prefix)_optim_best.h5", optim_best)
    save_result(
        "$(prefix)_tracked_candidates.h5",
        (
            candidates=Vector(tracked_candidates),
            fitnesses=Vector(tracked_fitnesses),
            threshold=fitness_threshold,
            n_tracked=length(tracked_candidates),
            trial_contribution=trial_contribution
        )
    )

    @info "Results saved to $(output_dir)"
end

"""
    load_calibration_results(output_dir, reef_id)

Load previously saved calibration results.
"""
function load_calibration_results(output_dir::String, reef_id::String)
    prefix = joinpath(output_dir, reef_id)

    # res = deserialize("$(prefix)_optim_state.dat")
    # optim_best = load_result("$(prefix)_optim_best.h5")
    # tracked_data = load_result("$(prefix)_tracked_candidates.h5")
    tmp = deserialize("$(prefix)_ensemble_output.dat")
    res = tmp.ensemble_res
    optim_best = tmp.best_params
    tracked_data = (;
        candidates=tmp.tracked_candidates,
        fitnesses=tmp.tracked_fitnesses,
        trial_contribution=tmp.trial_contribution,
        n_tracked=length(tmp.tracked_fitnesses)
    )

    @info "Loaded $(tracked_data.n_tracked) tracked candidates"

    return res,
    optim_best, tracked_data.candidates, tracked_data.fitnesses,
    tracked_data.trial_contribution
end

"""
    calculate_performance_metrics(sim, obs)

Calculate and return performance metrics for calibration.
"""
function calculate_performance_metrics(sim::Vector{Float64}, obs::Vector{Float64})
    # Error metrics
    init_mae = mean(abs(sim[1] - obs[1]))
    end_mae = mean(abs(sim[end] - obs[end]))
    rmse = Kora.RMSE(sim, obs)

    # Correlation
    kendall = Kora.kendall(sim, obs)
    pearson = Kora.pearson(sim, obs)

    # Bias
    u = abs((mean(sim) / mean(obs)) - 1.0)
    β = u / (1 + u)

    # Variability ratio
    std_obs = std(obs)
    std_hat = std(sim)
    _v = abs((std_hat / std_obs) - 1.0)
    α = (_v / (1 + _v))

    return (
        init_mae=init_mae,
        end_mae=end_mae,
        rmse=rmse,
        kendall=kendall,
        pearson=pearson,
        bias=β,
        variability=α
    )
end

function plot_calibration_results(
    reef_state::ReefState,
    env_conditions::AbstractDimArray,
    reef_obs::DataFrame,
    c_sim_indices::Vector{Int},
    c_ref_indices::Vector{Int},
    sim_year_range,
    cover::Vector{Float64},
    area::Float32,
    ecorrap_obs::DataFrame,
    output_file::String;
    v_sim_indices=nothing,
    v_ref_indices=nothing,
    ensemble_res=nothing,
    best_member_idx=nothing
)
    f = Figure(; size=(1400, 1000))  # Increased height for new subplot

    n_grps = n_groups(reef_state)
    n_ts = n_timesteps(reef_state)

    # Group cover for the best single member.
    # Prefer ensemble data when available (best_member_idx set); fall back to
    # reef_state (which may reflect an arbitrary run) otherwise.
    group_cover_best = if !isnothing(ensemble_res) && !isnothing(best_member_idx)
        (ensemble_res.group_cover[:, 1, :, best_member_idx] ./ area) .* 100.0
    else
        tmp = zeros(Float32, n_ts, n_grps)
        for ts in 1:n_ts, grp in 1:n_grps
            pop = coral_population(reef_state, ts, 1, grp)
            tmp[ts, grp] = sum(cover_cm_to_m2.(pop))
        end
        (tmp / area) * 100.0
    end

    # Very messy, to be cleaned up later.
    ecorrap_overlap_idx = indexin(ecorrap_obs.year, year.(sim_year_range))
    v_ecorrap = ecorrap_overlap_idx .∉ Ref(c_sim_indices)
    c_ecorrap = .!v_ecorrap

    # Get log pearson score of benthic observations (for one year in this case)
    v_bf_r_log = []
    v_ecorrap_obs = ecorrap_obs[v_ecorrap, :]
    for (i, r) in enumerate(eachrow(v_ecorrap_obs))
        _r = collect(values(r[2:end]))

        # Add small constant to handle any potential zero values
        v = vec(log.(group_cover_best[ecorrap_overlap_idx[v_ecorrap][i], :]')) .+ 1e-06
        push!(v_bf_r_log, cor(log.(_r), v))
    end
    v_bf_pearson = round(mean(v_bf_r_log); digits=2)

    c_bf_obs = reef_obs.MEAN_LIVE_CORAL[c_ref_indices]
    c_bf_rmse = round(
        Kora.RMSE(cover[c_sim_indices], c_bf_obs); digits=2
    )

    # Compute validation RMSE for best fit (post-disturbance manta tow years)
    v_bf_val_rmse = nothing
    if !isnothing(v_sim_indices) && !isempty(v_sim_indices)
        v_bf_obs_val = reef_obs.MEAN_LIVE_CORAL[v_ref_indices]
        v_bf_val_rmse = round(Kora.RMSE(cover[v_sim_indices], v_bf_obs_val); digits=2)
    end

    # Compute ensemble metrics and build annotation string
    sim_timeframe = Dates.year.(Date.(sim_year_range))
    if !isnothing(ensemble_res)
        ensemble_cover = (ensemble_res.cover[:, 1, :] / area) * 100.0
        ensemble_mean = vec(mean(ensemble_cover; dims=2))
        c_obs = reef_obs.MEAN_LIVE_CORAL[c_ref_indices]

        c_rmse = round(
            Kora.RMSE(ensemble_mean[c_sim_indices], c_obs); digits=2
        )
        c_pearson = round(
            Kora.pearson(ensemble_mean[c_sim_indices], c_obs); digits=2
        )

        v_aligned = ecorrap_overlap_idx[v_ecorrap]
        v_ecorrap_obs = ecorrap_obs[v_ecorrap, :]

        v_subset = ensemble_res.group_cover[v_aligned, 1, :, :] .+ 1e-06

        v_r_log = []
        for (i, r) in enumerate(eachrow(v_ecorrap_obs))
            _r = collect(values(r[2:end]))

            for a in axes(v_subset, 3)
                push!(v_r_log, cor(log.(_r), log.(v_subset[i, :, a])))
            end
        end

        mean_r_log = round(mean(v_r_log); digits=2)

        # Ensemble validation RMSE (post-disturbance manta tow years)
        v_ens_val_rmse = nothing
        if !isnothing(v_sim_indices) && !isempty(v_sim_indices)
            v_ens_obs_val = reef_obs.MEAN_LIVE_CORAL[v_ref_indices]
            v_ens_val_rmse = round(
                Kora.RMSE(ensemble_mean[v_sim_indices], v_ens_obs_val); digits=2
            )
        end

        bf_val_str = isnothing(v_bf_val_rmse) ? "" : " | Val RMSE: $(v_bf_val_rmse)"
        ens_val_str = isnothing(v_ens_val_rmse) ? "" : " | Val RMSE: $(v_ens_val_rmse)"
        metrics_annotation = L"""
            Best Single Member - Cal RMSE: %$(c_bf_rmse)%$(bf_val_str) | $r_{\mathrm{log}}$: %$(v_bf_pearson)\\
            Ensemble Mean - Cal RMSE: %$(c_rmse)%$(ens_val_str) | $r_{\mathrm{log}}$: %$(mean_r_log)
            """
    else
        bf_val_str = isnothing(v_bf_val_rmse) ? "" : " | Val RMSE: $(v_bf_val_rmse)"
        metrics_annotation = "Best Single Member - Cal RMSE: $(c_bf_rmse)$(bf_val_str) | r_log: $(v_bf_pearson)"
    end

    # Timeseries panel (total cover)
    ax1 = Axis(
        f[1, 1];
        xlabel="Date",
        ylabel="Total Coral Cover [%]",
        title="Ensemble Results",
        titlesize=20,
        xlabelsize=16,
        ylabelsize=16
    )

    # Shading regions for non-calibration years.
    # ecoRRAP validation indices (already computed above as ecorrap_overlap_idx[v_ecorrap])
    ecorrap_val_sim_idx = ecorrap_overlap_idx[v_ecorrap]

    # All simulation years between the first excluded post-disturbance year and the
    # last calibration year (exclusive) that are not themselves calibration years.
    # This includes years with no observation at all (e.g. 2019), which won't appear
    # in v_sim_indices but are still post-disturbance and should be shaded.
    post_dist_shade_years = if !isnothing(v_sim_indices) && !isempty(v_sim_indices)
        first_excl = minimum(sim_timeframe[v_sim_indices])
        last_calib = maximum(sim_timeframe[c_sim_indices])
        calib_years = Set(sim_timeframe[c_sim_indices])
        ecorrap_years = Set(sim_timeframe[ecorrap_val_sim_idx])
        [
            yr for yr in sim_timeframe
            if yr >= first_excl && yr < last_calib &&
               yr ∉ calib_years && yr ∉ ecorrap_years
        ]
    else
        Int[]
    end

    # ecoRRAP years — consecutive, so a single band; darker shade to distinguish
    ecorrap_shade_years = sim_timeframe[ecorrap_val_sim_idx]

    function shade_postdist_regions!(ax)
        for yr in post_dist_shade_years
            vspan!(ax, yr - 0.5, yr + 0.5; color=(:gray, 0.08))
        end
        if !isempty(ecorrap_shade_years)
            vspan!(
                ax,
                minimum(ecorrap_shade_years) - 0.5,
                maximum(ecorrap_shade_years) + 0.5;
                color=(:gray, 0.18)
            )
        end
    end

    shade_postdist_regions!(ax1)

    # Collect legend elements and labels
    legend_elements = []
    legend_labels = String[]

    # Plot ensemble trajectories if available (in background)
    if !isnothing(ensemble_res)
        # Ensemble member trajectories
        series!(
            ax1,
            sim_timeframe,
            ensemble_cover';
            solid_color=(:gray, 0.1)
        )

        # Create custom legend element for ensemble members
        ensemble_elem = LineElement(; color=(:gray, 0.3), linewidth=2)
        push!(legend_elements, ensemble_elem)
        push!(legend_labels, "Ensemble Members (n=$(size(ensemble_cover, 2)))")

        # Ensemble mean
        p_ens_mean = lines!(
            ax1,
            sim_timeframe,
            ensemble_mean;
            color=(:green, 0.9),
            linewidth=4
        )
        push!(legend_elements, p_ens_mean)
        push!(legend_labels, "Ensemble Mean")

        scatter!(
            ax1, sim_timeframe[c_sim_indices], ensemble_mean[c_sim_indices];
            color=(:orange, 0.8),
            markersize=8
        )
    end

    # Observations
    # draw band first so other elements are overlaid
    band!(
        ax1,
        reef_obs.SAMPLE_DATE,
        reef_obs.LOWER,
        reef_obs.UPPER;
        color=(:red, 0.1)
    )
    p_obs = scatterlines!(
        reef_obs.SAMPLE_DATE,
        reef_obs.MEAN_LIVE_CORAL;
        linestyle=:dashdotdot,
        color=(:red, 0.4),
        linewidth=2
    )
    push!(legend_elements, p_obs)
    push!(legend_labels, "Observations")

    scatter!(
        ax1, reef_obs.SAMPLE_DATE[c_ref_indices],
        reef_obs.MEAN_LIVE_CORAL[c_ref_indices];
        color=(:orange, 0.8),
        markersize=8
    )

    # Best single member trajectory
    p_best = scatterlines!(
        ax1, sim_timeframe, cover;
        linewidth=2,
        color=(:blue, 0.4)
    )
    push!(legend_elements, p_best)
    push!(legend_labels, "Best Single Member")

    p_matched = scatter!(
        ax1, sim_timeframe[c_sim_indices], cover[c_sim_indices];
        color=(:orange, 0.8),
        markersize=8
    )
    push!(legend_elements, p_matched)
    push!(legend_labels, "Matched Timepoints")

    ylims!(0.0, maximum(cover) + 10.0)

    # Add legend to the right
    Legend(f[1, 2], legend_elements, legend_labels; labelsize=16)

    # Metrics annotation inside the plot (bottom left)
    text!(ax1, 0.02, 0.03;
        text=metrics_annotation,
        align=(:left, :bottom),
        space=:relative,
        fontsize=16
    )

    hidexdecorations!(ax1; ticks=true, ticklabels=true, grid=false)

    # Group trajectories panel
    ax2 = Axis(
        f[2, 1];
        xlabel="Date",
        ylabel="Coral Cover by Group [%]",
        titlesize=20,
        xlabelsize=16,
        ylabelsize=16
    )
    shade_postdist_regions!(ax2)

    # Collect group legend elements and labels
    group_legend_elements = []
    group_legend_labels = String[]

    # Okabe-Ito colorblind-safe palette, paired with distinct linestyles/markers so
    # groups remain distinguishable in grayscale and under CVD simulation.
    FGROUP_COLOR = [
        "#E69F00",  # orange       - Tabular Acropora
        "#D55E00",  # vermillion   - Corymbose Acropora
        "#56B4E9",  # sky blue     - branching non-Acropora
        "#009E73",  # green        - Small massives
        "#0072B2",  # blue         - Large massives
    ]
    FGROUP_LINESTYLE = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    FGROUP_MARKER = [:circle, :rect, :utriangle, :diamond, :xcross]
    FLABELS = [
        "Tabular Acropora", "Corymbose Acropora",
        "branching non-Acropora", "Small massives", "Large massives"
    ]

    # Plot group trajectories with confidence intervals if ensemble available
    if !isnothing(ensemble_res)
        # Convert ensemble group cover to percentages
        ensemble_group_cover = (ensemble_res.group_cover[:, 1, :, :] / area) * 100.0

        for grp in n_grps:-1:1  # Reverse order for better layering
            # Calculate mean and 95% CI
            grp_mean = vec(mean(ensemble_group_cover[:, grp, :]; dims=2))
            grp_lower = [
                quantile(ensemble_group_cover[t, grp, :], 0.025) for t in 1:n_ts
            ]
            grp_upper = [
                quantile(ensemble_group_cover[t, grp, :], 0.975) for t in 1:n_ts
            ]

            # Plot confidence band
            band!(
                ax2,
                sim_timeframe,
                grp_lower,
                grp_upper;
                color=(FGROUP_COLOR[grp], 0.3)
            )

            # Plot mean line with distinct linestyle/marker for colorblind accessibility
            p_grp = scatterlines!(
                ax2,
                sim_timeframe,
                grp_mean;
                color=FGROUP_COLOR[grp],
                linestyle=FGROUP_LINESTYLE[grp],
                linewidth=2.5,
                marker=FGROUP_MARKER[grp],
                markersize=14,
                markercolor=FGROUP_COLOR[grp],
                strokewidth=1,
                strokecolor=:white
            )

            # `band!` does not support datetimes at this stage so we hide the xtick
            # labels to compensate.
            hidexdecorations!(ax2; ticks=true, ticklabels=true, grid=false)

            push!(group_legend_elements, p_grp)
            push!(group_legend_labels, FLABELS[grp])
        end
    else
        # Plot only best fit lines (no ensemble)
        for grp in n_grps:-1:1  # Reverse order for better layering
            p_grp = scatterlines!(
                ax2,
                sim_timeframe,
                group_cover_best[:, grp];
                color=FGROUP_COLOR[grp],
                linestyle=FGROUP_LINESTYLE[grp],
                linewidth=2.5,
                marker=FGROUP_MARKER[grp],
                markersize=14,
                markercolor=FGROUP_COLOR[grp],
                strokewidth=1,
                strokecolor=:white
            )
            push!(group_legend_elements, p_grp)
            push!(group_legend_labels, FLABELS[grp])
        end
    end

    # Add group legend to the right
    Legend(
        f[2, 2], reverse(group_legend_elements), reverse(group_legend_labels); labelsize=16
    )

    # DHW conditions - plot with dates instead of timesteps
    ax3 = Axis(
        f[3, 1];
        xlabel="Date",
        ylabel="DHW [°C-weeks]",
        titlesize=20,
        xlabelsize=16,
        ylabelsize=16
    )
    shade_postdist_regions!(ax3)

    # Get DHW data and create date vector
    dhw_data = env_conditions[:, :, At(:dhw)].data
    n_timesteps_dhw = size(dhw_data, 1)

    # Create date vector matching the DHW timesteps
    sim_year_range = Date.(sim_year_range)
    if length(sim_year_range) == n_timesteps_dhw
        dhw_dates = Dates.year.(sim_year_range)
    else
        # Fallback: create dates based on DHW data length
        start_date = first(sim_year_range)
        dhw_dates = start_date:Year(1):(start_date + Year(n_timesteps_dhw - 1))
    end

    # Plot DHW for each location
    for loc in axes(dhw_data, 2)
        lines!(ax3, dhw_dates, dhw_data[:, loc]; color=(:blue, 0.6))
    end

    # Link x-axes so zooming/panning is synchronized
    linkxaxes!(ax1, ax2, ax3)

    save(output_file, f; px_per_unit=DPI)
    @info "Plot saved to $(output_file)"

    return f
end

# ── Recruitment temporal sensitivity figures ──────────────────────────────────

"""Integer year ticks: every year for short series, otherwise Makie's default spacing."""
year_ticks(years) = length(years) <= 10 ? collect(years) : Makie.automatic

"""
Mark each disturbance year with a dashed red line, labelled at `ytop` when `label` is set.
`label_side` places the label to the `:right` or `:left` of the line.
"""
function mark_disturbances!(
    ax, disturbance_years, label::String, ytop; label_color=(:red, 0.8), label_side=:right
)
    halign, dx = label_side == :right ? (:left, 4) : (:right, -4)
    for yr in disturbance_years
        vlines!(ax, yr; color=(:red, 0.6), linewidth=1.5, linestyle=:dash)
        isempty(label) && continue
        text!(
            ax, yr, ytop; text=label, align=(halign, :top), offset=(dx, -2), fontsize=12,
            color=label_color
        )
    end
end

"""
    plot_recruitment_trajectories(years, with_recruitment, without_recruitment;
                                  band_quantiles=(0.25, 0.75), disturbance_years=Float64[],
                                  disturbance_label="", title="")

Two panels: (A) total coral cover with and without recruitment, (B) the difference
(recruitment-attributable cover). Lines are the mean across samples and bands span
`band_quantiles` (interquartile range by default). Inputs are `n_samples × n_years` cover
proportions (as returned by `recruitment_paired_runs`).
"""
function plot_recruitment_trajectories(
    years::AbstractVector,
    with_recruitment::AbstractMatrix,
    without_recruitment::AbstractMatrix;
    band_quantiles::Tuple{Float64,Float64}=(0.25, 0.75),
    disturbance_years=Float64[],
    disturbance_label::String="",
    title::String=""
)
    pct(Y) = 100 .* Y
    band(Y, p) = [quantile(pct(Y[:, t]), p) for t in axes(Y, 2)]
    mean_line(Y) = vec(mean(pct(Y); dims=1))
    q_lo, q_hi = band_quantiles

    colors = Makie.wong_colors()
    f = Figure(; size=(1400, 520))
    Label(f[0, 1:2], title; fontsize=20, font=:bold, tellwidth=false)

    ax_cover = Axis(f[1, 1]; xlabel="Year", ylabel="Coral cover (%)", title="(A) Total cover", xticks=year_ticks(years))
    for (Y, label, col) in [
        (with_recruitment, "With recruitment", colors[1]),
        (without_recruitment, "Without recruitment", colors[2])
    ]
        band!(ax_cover, years, band(Y, q_lo), band(Y, q_hi); color=(col, 0.2))
        lines!(ax_cover, years, mean_line(Y); color=col, linewidth=2.5, label=label)
    end
    axislegend(ax_cover; position=:lt, framevisible=false)

    diff = with_recruitment .- without_recruitment
    ax_diff = Axis(
        f[1, 2]; xlabel="Year", ylabel="Cover difference (percentage points)",
        title="(B) Recruitment-attributable cover (with − without)",
        xticks=year_ticks(years)
    )
    band!(ax_diff, years, band(diff, q_lo), band(diff, q_hi); color=(colors[3], 0.2))
    lines!(ax_diff, years, mean_line(diff); color=colors[3], linewidth=2.5)

    ytop_cover = maximum(max(maximum(band(Y, q_hi)), maximum(mean_line(Y)))
        for Y in (with_recruitment, without_recruitment))
    ytop_diff = max(maximum(band(diff, q_hi)), maximum(mean_line(diff)))
    mark_disturbances!(ax_cover, disturbance_years, disturbance_label, ytop_cover)
    mark_disturbances!(ax_diff, disturbance_years, disturbance_label, ytop_diff)

    return f
end

"""
    plot_recruitment_temporal_pawn(pawn_cube; scenarios=[:with, :without, :difference],
                                   disturbance_years=Float64[], disturbance_label="", title="")

Year × parameter heatmaps of dummy-normalized PAWN for each of `scenarios` in `pawn_cube`
(`factors × year × scenario`, scenarios `:with`, `:without`, `:difference`), one panel per
scenario. Rows are ordered by mean PAWN in the first scenario shown, and all panels share one
color scale; values at or below zero (not above the dummy) share the lowest color. Panels are
lettered only when more than one is shown (e.g. `scenarios=[:difference]` for a standalone
figure). `disturbance_label`, if set, is shown in the first panel only.
"""
function plot_recruitment_temporal_pawn(
    pawn_cube::AbstractDimArray;
    scenarios::Vector{Symbol}=[:with, :without, :difference],
    disturbance_years=Float64[],
    disturbance_label::String="",
    title::String=""
)
    factor_names = string.(collect(dims(pawn_cube, :factors)))
    years = collect(dims(pawn_cube, :year))
    n_factors = length(factor_names)
    n_panels = length(scenarios)

    shown = [parent(pawn_cube[scenario=At(s)]) for s in scenarios]
    row_order = sortperm(vec(mean(first(shown); dims=2)); rev=true)
    max_val = max(0.1, maximum(maximum, shown))

    panel_names = Dict(
        :with => "With recruitment",
        :without => "Without recruitment",
        :difference => "Recruitment-attributable (with − without)"
    )

    f = Figure(; size=(150 + 580 * n_panels, max(450, 32 * n_factors)))
    Label(f[0, 1:n_panels], title; fontsize=20, font=:bold, tellwidth=false)

    hm = nothing
    for (col, (scenario, vals)) in enumerate(zip(scenarios, shown))
        panel_title = n_panels > 1 ? "($('A' + col - 1)) $(panel_names[scenario])" :
            panel_names[scenario]
        ax = Axis(
            f[1, col];
            xlabel="Year",
            title=panel_title,
            xticks=year_ticks(years),
            yticks=(1:n_factors, factor_names[row_order]),
            yticklabelsvisible=col == 1
        )
        hm = heatmap!(
            ax, years, 1:n_factors, vals[row_order, :]';
            colormap=:viridis, colorrange=(0.0, max_val)
        )
        mark_disturbances!(
            ax, disturbance_years, col == 1 ? disturbance_label : "", n_factors + 0.5;
            label_color=:white, label_side=:left
        )
    end
    Colorbar(f[1, n_panels + 1], hm; label="PAWN index (relative to dummy)", width=14)

    return f
end
