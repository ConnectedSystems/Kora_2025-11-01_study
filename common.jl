using Revise, Infiltrator

using Serialization
using Random

using LinearAlgebra
using Distributions, Statistics, StatsBase
using KernelDensity
using CategoricalArrays

using Dates
using CSV, DataFrames
using NetCDF, YAXArrays
using Parquet2
import GeoDataFrames as GDF

using CairoMakie
using CoralFlow

Makie.inline!(true)

OUTPUT_DIR = "data"
FIG_DIR = "figs"
DPI = 300 / 96  # desired unit / pixels per inch

# Unnecessary/correlated factors to remove
surv_ignore_cols = [
    :class_train, :class_test, :class_train_mean, :class_test_mean, :surv_logclass,
    :logdiam, Symbol("days_t1.t2"), :cluster, :bleaching_scores,
    :surv, :diam, :survival_use, :growth_use, :class_test_std,
    :class_train_std, :depth_category, :dataset, :water_clarity,
    :site_new, :transition, :plot, :size, :sizenext,
    :clarified_note_2023_july, Symbol("clarified_note_2023.1"),
    :clarified_note_2021, :clarified_note_2022, :clarified_note_2023,
    :date_2021, :date_2022, :date_2023,
    :coral_cover_2021, :coral_cover_2022, :coral_cover_2023,
    :est_1yo_growth, :growth_rate
]

growth_ignore_cols = [
    :class_train, :class_test, :class_train_mean, :class_test_mean, :surv_logclass,
    :logdiam, Symbol("days_t1.t2"), :cluster, :bleaching_scores,
    :surv, :diamnext, :survival_use, :growth_use, :class_test_std,
    :class_train_std, :depth_category, :dataset, :water_clarity,
    :site_new, :transition, :plot, :logdiam, :growth, :lin_ext, :size, :sizenext,
    :clarified_note_2023_july, Symbol("clarified_note_2023.1"),
    :clarified_note_2021, :clarified_note_2022, :clarified_note_2023,
    :date_2021, :date_2022, :date_2023,
    :coral_cover_2021, :coral_cover_2022, :coral_cover_2023,
    :est_1yo_growth, :growth_rate
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
    for n in names(X)
        if eltype(X[!, n]) <: Union{AbstractString,Missing}
            X[!, n] .= Int64.(categorical(X[!, n]).refs)
        end

        if eltype(X[!, n]) <: Union{Float64,Missing}
            X[!, n] .= Float64.(X[!, n])
        end

        if eltype(X[!, n]) <: Union{Int64,Missing}
            X[!, n] .= Int64.(X[!, n])
        end
    end
end

const _DISPLAY_RENAMES = Dict(
    :Cscape_group => :functional_group,
    :diam_mort => :diameter,
    :temp => :temperature,
    :depth_cont => :depth,
    :plot_uid => :plot,
    :site_uid => :site,
    :ubed90_median => :bottom_stress
)

function rename_for_display!(df::DataFrame)
    cols = propertynames(df)
    pairs = [old => new for (old, new) in _DISPLAY_RENAMES if old in cols]
    isempty(pairs) || rename!(df, pairs...)
    return df
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

    cols = Dict{Symbol, Vector}(:Group => collect(fits.names))
    for m in metrics
        cols[Symbol("Train_$(m)")] = collect(Float64.(perf.train[m]))
        cols[Symbol("Test_$(m)")]  = collect(Float64.(perf.test[m]))
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
    CSV.write(joinpath(output_dir, "$(prefix)_coefficients.csv"), coeff_df)
end

include("src/sensitivity.jl")
include("src/parameter_assessment.jl")
include("src/calibration_structs.jl")
include("src/calibration_helpers.jl")