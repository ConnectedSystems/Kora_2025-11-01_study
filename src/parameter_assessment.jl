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
