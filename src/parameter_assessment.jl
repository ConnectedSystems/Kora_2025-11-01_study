using StatsBase

"""
Determine diversity of ensemble.

Checks how wide of an exploration the ensemble represents.

Mean distance:
    Average distance between all pairs of ensemble members.
    Higher values indicate greater diversity, lower values indicate parameter values are
    clustered.

Normalized Spread:
    Higher values indicate diverse exploration of parameter space, but equally, that
    equifinality is high (many solutions explain behavior)

    Low values indicate parameters are well-constrained or that exploration was not
    sufficiently wide.

Total variance:
    High values indicate variability in parameters
    Low values indicate the ensemble clusters in small region of parameter space.
"""
function effective_ensemble_diversity(ensemble_params)
    n = size(ensemble_params, 2)
    distances = [norm(ensemble_params[:, i] - ensemble_params[:, j])
                 for i in 1:n, j in 1:n]

    # Higher = more diverse
    mean_distance = mean(distances[distances .> 0])
    normalized_spread = mean(distances) / sqrt(size(ensemble_params, 1))

    # Higher = more diverse
    total_var = sum(var(ensemble_params; dims=2))

    # Normalize by parameter space diagonal

    return (
        mean_distance=mean_distance,
        normalized_spread=normalized_spread,
        total_variance=total_var
    )
end

"""
Determine constraint of each parameter.

Calculates Coefficient of Variation (stdev / mean), range ratio, and effective samples.

High CV and range indicate equifinality (spread of potential parameterizations).

CV: Low values indicate values are clustered around the mean (well-constrained)
Range ratio: Values closer to 0 indicate tight clusters (well-constrained)
             This is indicative, but sensitive to outliers.

"""
function parameter_identifiability_metrics(ensemble_params, param_names)
    n_params = length(param_names)
    metrics = DataFrame(;
        parameter=param_names,
        CV=zeros(length(param_names)), # coefficient of variation
        range_ratio=zeros(length(param_names)),  # coefficient of range
        MAD=zeros(length(param_names)),  # median absolute deviation
        rMAD=zeros(length(param_names))  # relative median absolute deviation
    )

    for p in 1:n_params
        values = ensemble_params[p, :]
        cv = std(values) / mean(values)

        max_val, min_val = maximum(values), minimum(values)
        range_ratio = (
            (max_val - min_val) / (max_val + min_val)
        )

        med_val = median(values)
        mad = median(abs.(values .- med_val))
        rmad = mad / med_val

        metrics[p, 2:end] = cv, range_ratio, mad, rmad
    end

    return metrics
end

"""
Per-candidate nearest-neighbour diversity check for an ensemble.

Unlike `effective_ensemble_diversity`'s single mean pairwise distance, which a
handful of duplicate clusters can barely move, this asks, for every candidate,
how close its closest neighbour is — surfacing near-duplicate solutions (repeated
optimizer convergence to the same basin) that an aggregate mean can hide.

Each parameter is rescaled to its `param_bounds` range before computing Euclidean
distance, so wide-ranged parameters (e.g. size means) don't dominate distances
over tightly-bounded ones (e.g. recruitment).

`ensemble_params` must be in shape [params ⋅ candidates], as with
`parameter_identifiability_metrics`. `param_bounds` is an iterable of `(lo, hi)`
pairs, one per row of `ensemble_params`/`param_names`.

`near_dup_threshold` is a fraction of the normalized parameter-space diagonal
(`sqrt(n_params)`); candidates whose nearest-neighbour distance falls below it are
flagged as near-duplicates.

Returns a per-candidate DataFrame with the nearest-neighbour index/distance.
"""
function parameter_nearest_neighbour_diversity(
    ensemble_params, param_names, param_bounds; near_dup_threshold=0.05
)
    n_params, n_candidates = size(ensemble_params)
    length(param_names) == n_params ||
        throw(ArgumentError("param_names must have one entry per row of ensemble_params"))

    lo = Float64.(first.(param_bounds))
    hi = Float64.(last.(param_bounds))
    span = hi .- lo
    normalized = (ensemble_params .- lo) ./ span

    nn_dist = fill(Inf, n_candidates)
    nn_idx = zeros(Int, n_candidates)
    for i in 1:n_candidates
        for j in 1:n_candidates
            i == j && continue
            d = norm(view(normalized, :, i) .- view(normalized, :, j))
            if d < nn_dist[i]
                nn_dist[i] = d
                nn_idx[i] = j
            end
        end
    end

    rel_nn_dist = nn_dist ./ sqrt(n_params)

    return DataFrame(;
        candidate=1:n_candidates,
        nearest_neighbour=nn_idx,
        nn_distance=nn_dist,
        rel_nn_distance=rel_nn_dist,
        near_duplicate=rel_nn_dist .< near_dup_threshold
    )
end

"""
Determine high correlations to identify trade-offs between parameters.

ensemble_params must be in shape [params ⋅ values]
"""
function parameter_correlation_analysis(ensemble_params, param_names; corr_threshold=0.7)
    # Calculate correlation matrix
    cor_matrix = cor(ensemble_params')

    # Identify strong correlations
    strong_correlations = []
    for i in 1:(length(param_names) - 1)
        for j in (i + 1):length(param_names)
            if abs(cor_matrix[i, j]) > corr_threshold
                push!(
                    strong_correlations,
                    (
                        Param1=param_names[i],
                        Param2=param_names[j],
                        Correlation=cor_matrix[i, j]
                    )
                )
            end
        end
    end

    return DataFrame(strong_correlations)
end

"""
    standardized_pearson(sim, obs)

Z-normalized Pearson correlation.

Z-score normalization removes absolute scale and centers the data, so only relative patterns matter.
"""
function standardized_pearson(sim, obs)
    return cor(zscore(sim), zscore(obs))
end

function log_pearson(sim, obs)
    return cor(log.(sim), log.(obs))
end
