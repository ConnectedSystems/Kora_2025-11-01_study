using Revise, Infiltrator

using Serialization
using Random

using LinearAlgebra
using Distributions, Statistics, StatsBase
using KernelDensity
using CategoricalArrays

using Dates
using CSV, DataFrames
using NCDatasets
using DimensionalData
using Parquet2
import GeoDataFrames as GDF

using CairoMakie
using Kora

# save_result / load_result — HDF5 storage for cached analysis results
include(joinpath(@__DIR__, "..", "src", "result_io.jl"))

OUTPUT_DIR = joinpath(@__DIR__, "..", "data")
FIG_DIR = joinpath(@__DIR__, "..", "figs")
EXT_DATA_DIR = joinpath(@__DIR__, "..", "..", "data")
DPI = 300 / 96  # desired unit / pixels per inch

# Explicit include-lists for the PAWN feature matrix. Chosen over a deny-list
# because a deny-list silently admits any new pipeline column into the
# sensitivity analysis unless someone remembers to blacklist it — this nearly
# happened with `depth_gapfilled` itself. An include-list fails safe: a new
# column is simply absent from PAWN until someone deliberately adds it here.
#
# Verified via a live `Kora.process_ecorrap_models` run (offshore_north)
# against the current EcoRRAP data.
# `wave_hs_mean` is deliberately excluded here (sparse coverage) but is added
# conditionally by 01c's wave-subset analysis.
growth_include_cols = [
    :diam, :depth_gapfilled, :ereefs_temp_max, :ereefs_temp_mean,
    :functional_group, :habitat_type, :plot_lat, :plot_lon, :taxa, :wave_ubed90,
]

surv_include_cols = [
    :diam_mort, :depth_gapfilled, :ereefs_temp_max, :ereefs_temp_mean,
    :functional_group, :habitat_type, :plot_lat, :plot_lon, :taxa, :wave_ubed90,
]

ENSEMBLE_PARAM_NAMES = [
    "Density",                      # 1
    "Prop: Tab Acro",              # 2
    "Prop: Cor Acro",              # 3
    "Prop: Cor non-Acro",          # 4
    "Prop: Sm Mass",               # 5
    "Prop: Lg Mass",               # 6
    "Size μ: Tab Acro",            # 7
    "Size μ: Cor Acro",            # 8
    "Size μ: Cor non-Acro",        # 9
    "Size μ: Sm Mass",             # 10
    "Size μ: Lg Mass",             # 11
    "Size σ: Tab Acro",            # 12
    "Size σ: Cor Acro",            # 13
    "Size σ: Cor non-Acro",        # 14
    "Size σ: Sm Mass",             # 15
    "Size σ: Lg Mass",             # 16
    "Scaler: Tab Acro",            # 17
    "Scaler: Cor Acro",            # 18
    "Scaler: Cor non-Acro",        # 19
    "Scaler: Sm Mass",             # 20
    "Scaler: Lg Mass",             # 21
    "External Recruitment",        # 22
    "Self-seeding"                 # 23
]

function cleanup_features!(X::DataFrame)
    drop_cols = String[]
    for n in names(X)
        T = nonmissingtype(eltype(X[!, n]))
        if T <: AbstractString
            X[!, n] = Float64.(categorical(X[!, n]).refs)
        elseif T <: AbstractFloat
            # missing → NaN
            X[!, n] = [ismissing(x) ? NaN : Float64(x) for x in X[!, n]]
        elseif T <: Real  # Bool, Int32, Int64, etc.
            # missing → 0.0
            X[!, n] = [ismissing(x) ? 0.0 : Float64(x) for x in X[!, n]]
        else
            push!(drop_cols, n)
        end
    end
    if !isempty(drop_cols)
        @warn "Dropping non-numeric columns from feature matrix" drop_cols
        select!(X, Not(drop_cols))
    end
end

const _DISPLAY_RENAMES = Dict(
    :Cscape_group         => :functional_group,
    :cscape_group         => :functional_group,
    :diam                 => :diameter,
    :diam_mort            => :diameter,
    :depth_gapfilled      => Symbol("Gapfilled Depth (m)"),
    :plot_uid             => :plot,
    :reef_habitat         => :site,
    :ubed90_median        => Symbol("Bottom Stress"),
    :temp                 => :ipm_temperature,
    :temp_max_mean        => :logger_temperature_mean_max,
    :temp_mean_mean       => :logger_temperature_mean_mean,
    :habitat_type         => :habitat,
    # Environmental covariates — explicit human-readable labels
    :wave_ubed90          => Symbol("Bottom Stress"),
    :depth_m              => Symbol("Depth (m)"),
    :depth_cat            => Symbol("Depth Category"),
    :ereefs_temp_mean     => Symbol("eReefs Temperature Mean"),
    :ereefs_temp_max      => Symbol("eReefs Temperature Max"),
    # Identifier / classification columns
    :colony_id            => :colony,
    :lat                  => :latitude,
    :lon                  => :longitude,
    :plot_lat             => :latitude,
    :plot_lon             => :longitude,
    :taxa                 => :taxon,
)

function rename_for_display!(df::DataFrame)
    cols = propertynames(df)
    # Skip a pair if its target name is already taken by another column --
    # the source column is left as-is and picked up by the underscore
    # fallback pass below instead of erroring in rename!. Guards against
    # e.g. a future column colliding with `habitat_type => habitat`.
    occupied = Set(string.(cols))
    pairs = Pair{Symbol,Symbol}[]
    for (old, new) in _DISPLAY_RENAMES
        old in cols || continue
        string(new) ∈ occupied && continue
        push!(pairs, old => new)
        push!(occupied, string(new))
    end
    isempty(pairs) || rename!(df, pairs...)

    # Fallback: any column name still containing underscores gets underscores
    # replaced with spaces and titlecased, so no snake_case labels appear in
    # figures.  Skip if the candidate name already exists as another column.
    occupied = Set(string.(propertynames(df)))
    fallback_pairs = Pair{Symbol,Symbol}[]
    for col in propertynames(df)
        s = string(col)
        occursin('_', s) || continue
        candidate = Symbol(titlecase(replace(s, "_" => " ")))
        string(candidate) ∈ occupied && continue
        push!(fallback_pairs, col => candidate)
        push!(occupied, string(candidate))
    end
    isempty(fallback_pairs) || rename!(df, fallback_pairs...)

    return df
end

"""
    plot_pawn_heatmap(Si, title; fig_size)

Create a PAWN sensitivity heatmap from a `pawn()` result slice (already filtered to
desired stats via `[PAWNᵢ=At(...)]`), sorting factors by mean PAWN index so the
most influential factor appears at the top.

Axis labels are always set explicitly so factor names are never left to chance.
"""
function plot_pawn_heatmap(
    Si::AbstractDimArray,
    title::String;
    stats::Vector{Symbol}=[:mean, :std],
    fig_size::Tuple{Int,Int}=(900, 420),
    xticklabelrotation::Real=π / 4
)
    # Sort factors by mean PAWN index — slice scalar At() to avoid At(vector) ambiguity
    factor_order = sortperm(collect(Si[PAWNᵢ=At(:mean)]); rev=true)
    factor_labels = string.(collect(dims(Si, 1)))[factor_order]

    # Build data matrix by stacking individual stat slices (n_factors × n_stats)
    data = hcat([collect(Si[PAWNᵢ=At(s)])[factor_order] for s in stats]...)
    stat_labels = string.(stats)

    f = Figure(; size=fig_size)
    ax = Axis(f[1, 1])
    hm = heatmap!(ax, data; colorrange=(-0.1, max(maximum(data), 0.1)), colormap=:viridis)
    Colorbar(f[1, 2], hm; label="PAWN Index")

    ax.xticks = (1:length(factor_labels), factor_labels)
    ax.yticks = (1:length(stat_labels), stat_labels)
    ax.yreversed = true
    ax.xticklabelrotation = xticklabelrotation
    ax.xlabelsize = 14
    ax.ylabelsize = 14
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12
    ax.titlesize = 14
    ax.title = title

    resize_to_layout!(f)
    return f
end

"""
    plot_nn_diversity_histogram(nn_df, title; near_dup_threshold)

Histogram of per-candidate nearest-neighbour distances (from
`parameter_nearest_neighbour_diversity`), with `near_dup_threshold` marked so
clusters of near-duplicate ensemble members are visible at a glance.
"""
function plot_nn_diversity_histogram(
    nn_df::DataFrame,
    title::String;
    near_dup_threshold::Real=0.05,
    fig_size::Tuple{Int,Int}=(700, 400)
)
    f = Figure(; size=fig_size)
    ax = Axis(
        f[1, 1];
        xlabel="Nearest-neighbour distance (fraction of normalized param-space diagonal)",
        ylabel="Count",
        title=title
    )
    hist!(ax, nn_df.rel_nn_distance; bins=40, color=(:steelblue, 0.7))
    vlines!(ax, near_dup_threshold; color=:red, linestyle=:dash, linewidth=1.5)
    text!(
        ax, near_dup_threshold, 0.0;
        text="near-duplicate\nthreshold", color=:red, fontsize=11,
        align=(:left, :bottom), offset=(4, 4)
    )

    return f
end

"""
    export_model_summaries(fits, output_dir, prefix)

Save two CSV files to `output_dir`:
  - `<prefix>_performance.csv`  — train/test metrics per group
  - `<prefix>_coefficients.csv` — polynomial coefficients per group
"""
function export_model_summaries(fits, output_dir::String, prefix::String)
    mkpath(output_dir)

    # ── Performance ───────────────────────────────────────────────────────────
    perf = fits.performance
    metrics = keys(perf.train)  # e.g. (:RMSE, :R2, :pearson, :spearman, :kendall)

    cols = Dict{Symbol,Vector}(:Group => collect(fits.names))
    for m in metrics
        cols[Symbol("Train_$(m)")] = collect(Float64.(perf.train[m]))
        cols[Symbol("Test_$(m)")] = collect(Float64.(perf.test[m]))
    end

    # Interleave Train/Test columns for readability
    col_order = [:Group]
    for m in metrics
        push!(col_order, Symbol("Train_$(m)"), Symbol("Test_$(m)"))
    end

    perf_df = DataFrame(cols)[!, col_order]
    CSV.write(joinpath(output_dir, "$(prefix)_performance.csv"), perf_df)

    # ── Coefficients ──────────────────────────────────────────────────────────
    polys = getfield.(fits.models, :poly)
    coeff_df = DataFrame(;
        Group=collect(fits.names),
        Model=string.(polys)
    )
    return CSV.write(joinpath(output_dir, "$(prefix)_coefficients.csv"), coeff_df)
end

const _PARQUET_TYPES = Union{Real,AbstractString,Date,DateTime}

"""
Coerce columns Parquet2 cannot represent (and all-`missing` columns) to strings.
"""
function _parquet_safe(df::DataFrame)::DataFrame
    out = copy(df)
    for n in names(out)
        T = nonmissingtype(eltype(out[!, n]))
        if T === Union{} || T === Any || !(T <: _PARQUET_TYPES)
            out[!, n] = [ismissing(x) ? missing : string(x) for x in out[!, n]]
        end
    end

    return out
end

"""
    export_fit_data(groupings, output_dir, prefix)

Save the exact data subset the models were fitted to, as
`<prefix>_fitdata.parquet` and `<prefix>_fitdata.csv`.

`groupings` is the `OrderedDict{String, DataFrame}` returned by
`process_growth_models` / `process_survival_models`. Groups are stacked into a
single table with a leading `functional_group` column holding the grouping key,
so every row is traceable to the model it informed. The train/test split
columns (`class_train`, `class_test`, `class_*_mean`, ...) are retained, as the
fits cannot be reproduced without them.

Three columns would otherwise state the same group: `Cscape_group` (set to the
grouping key on every row by `collate_functional_groups`) and the dataset's own
`functional_group` (an incomplete upstream copy -- populated for ~92% of rows,
never in disagreement). Both are dropped in favour of the key itself, under the
preferred name.
"""
function export_fit_data(groupings, output_dir::String, prefix::String)
    mkpath(output_dir)

    frames = DataFrame[]
    for (grp, df) in groupings
        isempty(df) && continue
        sub = copy(df)

        redundant = intersect(
            [:functional_group, :Cscape_group, :cscape_group], propertynames(sub)
        )
        isempty(redundant) || select!(sub, Not(redundant))

        insertcols!(sub, 1, :functional_group => fill(grp, nrow(sub)))
        push!(frames, sub)
    end

    if isempty(frames)
        @warn "No data to export for $(prefix)"
        return nothing
    end

    combined = vcat(frames...; cols=:union)

    csv_path = joinpath(output_dir, "$(prefix)_fitdata.csv")
    CSV.write(csv_path, combined)

    Parquet2.writefile(
        joinpath(output_dir, "$(prefix)_fitdata.parquet"), _parquet_safe(combined)
    )

    return csv_path
end

include(joinpath(@__DIR__, "..", "src", "sensitivity.jl"))
include(joinpath(@__DIR__, "..", "src", "parameter_assessment.jl"))
include(joinpath(@__DIR__, "..", "src", "calibration_structs.jl"))
include(joinpath(@__DIR__, "..", "src", "calibration_helpers.jl"))