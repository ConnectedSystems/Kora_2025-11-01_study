"""
Produce annual mean benthic cover estimates by functional group for Moore Reef (ONMO)
and Masig Reef (TSMA) from the expanded EcoRRAP photoquadrat data.

Reads cover_estimate_DESCRIPTION.csv (one row per transect, species as columns) and
aggregates to annual mean cover per functional group by averaging across transects.

Outputs:
  data/ecorrap_benthic/moore_estimate.csv
  data/ecorrap_benthic/masig_estimate.csv
  data/Masig_Reef_EcoRRAP_estimate.csv   (total cover summary used by calibration scripts)
"""

using CSV
using DataFrames
using Statistics

include("./common.jl")

const BENTHIC_FILE = joinpath(
    OUTPUT_DIR, "ecorrap_benthic", "latest", "cover_estimate_DESCRIPTION.csv"
)
const LABELSET_FILE = joinpath(
    OUTPUT_DIR, "ecorrap_benthic", "EcoRRAP_Labelset_mapping_Maren_Toor_2026-01-20.csv"
)

const GROUP_ORDER = [
    "Tabular Acropora",
    "Corymbose Acropora",
    "branching non-Acropora",
    "Small massive",
    "Large massive"
]

# ─── Load data ────────────────────────────────────────────────────────────────

benthic = CSV.read(BENTHIC_FILE, DataFrame; missingstring=["", "NA", "N/A"])

labelset = CSV.read(LABELSET_FILE, DataFrame; missingstring="NA")
mapped = unique(labelset[:, [:DESCRIPTION, :MAPPING]])
mapped = mapped[.!ismissing.(mapped.MAPPING), :]

# ─── Build species → functional group lookup ──────────────────────────────────
# DESCRIPTION values in the labelset use dot notation ("Acropora.table").
# The benthic CSV column names may use spaces; build a map covering both forms.

meta_cols = Set([:survey_title, :year, :site, :depth_m, :site_reef_name, :transect])
cover_cols = [
    n for n in propertynames(benthic)
    if n ∉ meta_cols && !endswith(string(n), "_STE")
]

# Map each cover column to its functional group (if any)
col_to_group = Dict{Symbol,String}()
for r in eachrow(mapped)
    desc_dot = string(r.DESCRIPTION)
    desc_space = replace(desc_dot, "." => " ")
    for col in cover_cols
        s = string(col)
        if s == desc_dot || s == desc_space
            col_to_group[col] = string(r.MAPPING)
        end
    end
end

n_matched = length(col_to_group)
@info "Labelset matched $n_matched / $(length(cover_cols)) cover columns"

# ─── Estimate function ─────────────────────────────────────────────────────────

"""
Compute annual mean cover by functional group for a single site.
Each transect row is first reduced to per-group cover (sum of constituent species),
then averaged across transects within each year.
"""
function site_annual_cover(df::DataFrame, site_code::String)::DataFrame
    site_df = df[df.site .== site_code, :]
    isempty(site_df) && error("No rows found for site $site_code")

    rows = NamedTuple[]
    for grp in groupby(site_df, :year)
        year = grp.year[1]

        transect_covers = Dict(g => Float64[] for g in GROUP_ORDER)

        for row in eachrow(grp)
            for (col, fg) in col_to_group
                fg ∈ GROUP_ORDER || continue
                v = row[col]
                ismissing(v) && continue
                push!(transect_covers[fg], Float64(v))
            end
        end

        # Average across transects; functional group total per transect is the
        # sum of its constituent species, so we accumulate per-species values
        # and take the mean across transects.  Where a group has no data, emit 0.
        row_vals = Any[year]
        for g in GROUP_ORDER
            vals = transect_covers[g]
            push!(row_vals, isempty(vals) ? 0.0 : mean(vals))
        end

        push!(rows, NamedTuple{(:year, Symbol.(GROUP_ORDER)...)}(Tuple(row_vals)))
    end

    isempty(rows) && error("No data aggregated for site $site_code")
    result = DataFrame(rows)
    sort!(result, :year)
    return result
end

# ─── Moore Reef ───────────────────────────────────────────────────────────────

@info "Processing Moore Reef (ONMO)..."
moore_cover = site_annual_cover(benthic, "ONMO")

@info moore_cover
CSV.write(joinpath(OUTPUT_DIR, "ecorrap_benthic", "moore_estimate.csv"), moore_cover)
@info "Written: moore_estimate.csv"

# ─── Masig Reef ───────────────────────────────────────────────────────────────

@info "Processing Masig Reef (TSMA)..."
masig_cover = site_annual_cover(benthic, "TSMA")

@info masig_cover
CSV.write(joinpath(OUTPUT_DIR, "ecorrap_benthic", "masig_estimate.csv"), masig_cover)
@info "Written: masig_estimate.csv"

# Total cover summary consumed by calibration scripts (03b / 04a)
masig_total = DataFrame(;
    report_year=masig_cover.year,
    mean=map(sum, eachrow(masig_cover[:, GROUP_ORDER])),
    lower=map(sum, eachrow(masig_cover[:, GROUP_ORDER])),
    upper=map(sum, eachrow(masig_cover[:, GROUP_ORDER]))
)
CSV.write(joinpath(OUTPUT_DIR, "Masig_Reef_EcoRRAP_estimate.csv"), masig_total)
@info "Written: Masig_Reef_EcoRRAP_estimate.csv"
