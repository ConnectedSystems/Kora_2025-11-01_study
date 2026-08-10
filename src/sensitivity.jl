using StaticArrays
using Logging

using Random
using Statistics
using Distributions
using HypothesisTests
using HypothesisTests: ApproximateKSTest

import Distributions: sample
import QuasiMonteCarlo as QMC
import QuasiMonteCarlo: SobolSample, OwenScramble

using DimensionalData

"""
    DataCube(data::AbstractArray; kwargs...)::DimArray

Constructor for DimArray. When used with `axes_names`, the axes labels will be UnitRanges
from 1 up to that axis length.

# Arguments
- `data` : Array of data to be used when building the DimArray
- `axes_names` : Tuple of axes names
- `properties` : NamedTuple of metadata to be added to the DimArray
"""
function DataCube(
    data::AbstractArray; properties::Dict{Symbol,Any}=Dict{Symbol,Any}(), kwargs...
)::DimArray
    return DimArray(
        data, Tuple(Dim{name}(val) for (name, val) in kwargs); metadata=properties
    )
end
function DataCube(
    data::AbstractArray, axes_names::Tuple; properties::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::DimArray
    return DataCube(
        data; properties=properties, NamedTuple{axes_names}(1:len for len in size(data))...
    )
end
function DataCube(
    data::AbstractArray, axes_names::Tuple, properties::Dict{Symbol,Any}
)::DimArray
    return DataCube(
        data; properties=properties, NamedTuple{axes_names}(1:len for len in size(data))...
    )
end

"""
    ZeroDataCube(; T::DataType=Float64, kwargs...)::DimArray
    ZeroDataCube(axes_names::Tuple, axes_sizes::Tuple; T::DataType=Float64)::DimArray

Constructor for DimArray with all entries equal zero. When `axes_name` and `axes_sizes`
are passed, all axes labels will be ranges.

# Arguments
- `axes_names` : Tuple of axes names
- `axes_sizes` : Tuple of axes sizes
- `properties` : NamedTuple of metadata to be added to the DimArray
"""
function ZeroDataCube(;
    T::Type{D}=Float64, properties::Dict{Symbol,Any}=Dict{Symbol,Any}(), kwargs...
)::DimArray where {D}
    return DataCube(
        zeros(T, [length(val) for (name, val) in kwargs]...); properties=properties,
        kwargs...
    )
end
function ZeroDataCube(
    axes_names::Tuple,
    axes_sizes::Tuple;
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    T::Type{D}=Float64
)::DimArray where {D}
    return ZeroDataCube(;
        T=T, properties=properties, NamedTuple{axes_names}(1:size for size in axes_sizes)...
    )
end
function ZeroDataCube(
    axes_names::Tuple,
    axes_sizes::Tuple,
    properties::Dict{Symbol,Any};
    T::Type{D}=Float64
)::DimArray where {D}
    return ZeroDataCube(;
        T=T, properties=properties, NamedTuple{axes_names}(1:size for size in axes_sizes)...
    )
end

"""
    ks_statistic(ks)

Calculate the Kolmogorov-Smirnov test statistic.
"""
function ks_statistic(ks::ApproximateKSTest)::Float64
    n::Float64 = (ks.n_x * ks.n_y) / (ks.n_x + ks.n_y)

    return sqrt(n) * ks.δ
end

"""
    pawn(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, factor_names::Vector{String}; S::Int64=10)::DimArray
    pawn(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::DimArray
    pawn(X::AbstractDimArray, y::Union{AbstractDimArray,AbstractVector{<:Real}}; S::Int64=10)::DimArray
    pawn(X::Union{DataFrame,AbstractMatrix{<:Real}}, y::AbstractMatrix{<:Real}; S::Int64=10)::DimArray

Calculates the PAWN sensitivity index.
Implementation and docstring below adapted from the SALib Python package.

The PAWN method (by Pianosi and Wagener) is a moment-independent approach to Global
Sensitivity Analysis. Outputs are characterized by their Cumulative Distribution Function
(CDF), quantifying the variation in the output distribution after conditioning an input over
"slices" (\$S\$) - the conditioning intervals. If both distributions coincide at all slices
(i.e., the distributions are similar or identical), then the factor is deemed
non-influential.

This implementation applies the Kolmogorov-Smirnov test as the distance measure and returns
summary statistics (min, lower bound, mean, median, upper bound, max, std, and cv) over the
slices. The statistics are also z-score normalized using results for a dummy factor as the
minimum threshold (i.e., `(x - dummy) / stdev(x)`).

# Arguments
- `rs` : ResultSet
- `X` : Model inputs
- `y` : Model outputs
- `factor_names` : Names of each factor represented by columns in `X`
- `S` : Number of slides (default: 10)

# Returns
DimArray, of min, mean, lower bound, median, upper bound, max, std, and cv summary statistics.

# References
1. Pianosi, F., Wagener, T., 2018.
   Distribution-based sensitivity analysis from a generic input-output sample.
   Environmental Modelling & Software 108, 197-207.
   https://doi.org/10.1016/j.envsoft.2018.07.019

2. Baroni, G., Francke, T., 2020.
   GSA-cvd
   Combining variance- and distribution-based global sensitivity analysis
   https://github.com/baronig/GSA-cvd

3. Puy, A., Lo Piano, S., & Saltelli, A. 2020.
   A sensitivity analysis of the PAWN sensitivity index.
   Environmental Modelling & Software, 127, 104679.
   https://doi.org/10.1016/j.envsoft.2020.104679

4. https://github.com/SAFEtoolbox/Miscellaneous/blob/main/Review_of_Puy_2020.pdf

# Extended help
Pianosi and Wagener have made public their review responding to a critique of their method
by Puy et al., (2020). A key criticism by Puy et al. was that the PAWN method is sensitive
to its tuning parameters and thus may produce biased results. The tuning parameters referred
to are the number of samples (\$N\$) and the number of conditioning points - \$n\$ in Puy et
al., but denoted as \$S\$ here.

Puy et al., found that the ratio of \$N\$ (number of samples) to \$S\$ has to be
sufficiently high (\$N/S > 80\$) to avoid biased results. Pianosi and Wagener point out this
requirement is not particularly difficult to meet. Using the recommended value
(\$S := 10\$), a sample of 1024 runs (small for purposes of Global Sensitivity Analysis)
meets this requirement (\$1024/10 = 102.4\$). Additionally, lower values of \$N/S\$ is more
an indication of faulty experimental design moreso than any deficiency of the PAWN method.
"""
function pawn(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real},
    factor_names::Vector{String};
    S::Int64=10
)::DimArray
    N, D = size(X)
    step = 1 / S
    seq = 0.0:step:1.0

    X = copy(X)
    X = hcat(X, randn(N))

    # Add 1 dimension for dummy factor
    D = D + 1

    # Preallocate result structures
    X_q = @MVector zeros(S + 1)
    pawn_t = @MArray zeros(S, D)
    results = @MArray zeros(D, 8)
    q_stats = [0.025, 0.5, 0.975]

    # Hide warnings from HypothesisTests
    with_logger(NullLogger()) do
        for d_i in 1:D
            X_di = @view(X[:, d_i])
            X_q .= quantile(X_di, seq)

            Y_sel = @view(y[X_q[1] .<= X_di .<= X_q[2]])
            if length(Y_sel) > 0
                pawn_t[1, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            for s in 2:S
                Y_sel = @view(y[X_q[s] .< X_di .<= X_q[s + 1]])
                if length(Y_sel) == 0
                    continue  # no available samples
                end

                pawn_t[s, d_i] = ks_statistic(ApproximateTwoSampleKSTest(Y_sel, y))
            end

            p_ind = @view(pawn_t[:, d_i])
            p_mean = mean(p_ind)
            p_sdv = std(p_ind)
            p_cv = p_sdv ./ p_mean
            p_lb, p_med, p_ub = quantile(p_ind, q_stats)
            results[d_i, :] .= (
                minimum(p_ind),
                p_lb,
                p_mean,
                p_med,
                p_ub,
                maximum(p_ind),
                p_sdv,
                p_cv
            )
        end
    end

    replace!(results, NaN => 0.0, Inf => 0.0)

    # Range normalize, relative to dummy factor such that negative values indicate
    # factor was below insensitive threshold as indicated by the dummy factor.
    dummy = results[end, :]
    tmp = @view results[1:(end - 1), :]
    results = (tmp .- dummy') ./ std(tmp)

    # Handle zero division (because tmp - dummy can resolve to zero)
    replace!(results, NaN => 0.0, Inf => 0.0)

    col_names = [:min, :lb, :mean, :median, :ub, :max, :std, :cv]
    row_names = Symbol.(factor_names)
    return DataCube(results; factors=row_names, PAWNᵢ=col_names)
end
function pawn(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::DimArray
    return pawn(Matrix(X), y, names(X); S=S)
end
function pawn(
    X::AbstractDimArray,
    y::AbstractVector{<:Real};
    S::Int64=10
)::DimArray
    return pawn(parent(X), y, string.(collect(dims(X, 2))); S=S)
end
function pawn(
    X::AbstractDimArray,
    y::AbstractDimArray;
    S::Int64=10
)::DimArray
    # Boolean indexing in pawn errors when the mask selects nothing, so vec(y) is required
    return pawn(parent(X), vec(y), string.(collect(dims(X, 2))); S=S)
end
function pawn(
    X::Union{DataFrame,AbstractDimArray},
    y::AbstractMatrix{<:Real};
    S::Int64=10
)::DimArray
    N, D = size(y)
    if N > 1 && D > 1
        msg::String = string(
            "The current implementation of PAWN can only assess a single quantity",
            " of interest at a time."
        )
        throw(ArgumentError(msg))
    end

    # The wrapped call to `vec()` handles cases where matrix-like data type is passed in
    # (N x 1 or 1 x D) and so ensures a vector is passed along
    return pawn(X, vec(y); S=S)
end
