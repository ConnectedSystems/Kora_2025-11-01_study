"""
Collate the expanded EcoRRAP dataset by joining:
  1. Individual coral tracking data (IPM)
  2. Community benthic cover aggregated to functional groups
  3. Oceanographic period statistics (output of 00a_prep_oceanographic_stats.jl)

Join keys:
  benthic:       (site_code, habitat_area, depth_cat, survey_year)
  oceanographic: (ocn_site_code, depth_cat, survey_year)

survey_year convention: May–April, labelled by May start year.
  Transition "2021_2022" → survey_year 2021 (surveys occur ~March–April, before May).

Output: data/ecorrap_expanded.parquet
"""

using CSV
using DataFrames
using Parquet2
using Statistics
using Dates
using Kora

# ─── Paths ────────────────────────────────────────────────────────────────────

const DATA_DIR = joinpath(@__DIR__, "..", "data")

const IPM_FILE = joinpath(DATA_DIR, "ecorrap_adult_juv_combined_2021_2023_24062025.csv")
const JUV_FILE = joinpath(DATA_DIR, "ecorrap_juv_data_2021_2023_24062025.csv")
const BENTHIC_FILE = joinpath(
    DATA_DIR, "ecorrap_benthic", "latest", "cover_estimate_DESCRIPTION.csv"
)
const OCN_STATS_FILE = joinpath(DATA_DIR, "ecorrap_oceanographic", "ocn_annual_stats.parquet")
const SPECIES_FILE = joinpath(DATA_DIR, "ecorrap_to_cscape_species.csv")
const LABELSET_FILE = joinpath(
    DATA_DIR, "ecorrap_benthic",
    "EcoRRAP_Labelset_mapping_Maren_Toor_2026-01-20.csv"
)
const LOOKUP_FILE = joinpath(DATA_DIR, "reef_site_code_lookup.csv")
const OUTPUT_FILE = joinpath(DATA_DIR, "ecorrap_expanded.parquet")

# ─── Constants ────────────────────────────────────────────────────────────────

const FUNCTIONAL_GROUPS = [
    "acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive"
]

# Maps the MAPPING display names in the labelset CSV to internal functional group names
const LABELSET_MAPPING_TO_FG = Dict(
    "Tabular Acropora" => "acro_table",
    "Corymbose Acropora" => "acro_corym",
    "branching non-Acropora" => "corym_non_acro",
    "Small massive" => "small_massive",
    "Large massive" => "large_massive"
)

# Maps the IPM SITE field to the short habitat-area code used in benthic survey_title
const SITE_TO_HABITAT = Dict(
    "Back1" => "BA1", "Back2" => "BA2", "Back3" => "BA3",
    "Back4" => "BA4", "Back5" => "BA5", "Back6" => "BA6",
    "Front1" => "FR1", "Front2" => "FR2",
    "Flank1" => "FL1", "Flank2" => "FL2",
    "Lagoon1" => "LA1", "Lagoon2" => "LA2"
)

# ─── Helpers ──────────────────────────────────────────────────────────────────

"""Return the survey_year from a transition string like "2021_2022" → 2021."""
function transition_to_survey_year(t::AbstractString)::Union{Int,Missing}
    parts = split(t, "_")
    length(parts) < 1 && return missing
    v = tryparse(Int, parts[1])
    return isnothing(v) ? missing : v
end

"""
Parse (habitat_area, depth_cat) from a survey_title field.
Format: SITE_HABITATD_P#_YEAR_DEPTHMETRE  e.g. ONMO_BA1D_P1_2021_4M
The depth indicator (D/S) is the last character of the second underscore-delimited token.
Returns ("", "") on parse failure.
"""
function parse_survey_title(title::AbstractString)::Tuple{String,String}
    parts = split(title, "_")
    length(parts) < 2 && return ("", "")
    hab = string(parts[2])
    length(hab) < 2 && return (hab, "")
    return (hab[1:(end - 1)], string(last(hab)))   # ("BA1", "D")
end

# ─── Species → functional group maps ─────────────────────────────────────────

"""
Build a Dict mapping DESCRIPTION column name → internal functional group name,
from the EcoRRAP labelset CSV (DESCRIPTION → MAPPING, one-to-one, NA rows skipped).
DESCRIPTION values use dot notation ("Acropora.table"); the map includes both the
raw dot form and the space form so either column naming convention is matched.
"""
function build_benthic_fg_map(labelset_df::DataFrame)::Dict{String,String}
    result = Dict{String,String}()
    for r in eachrow(labelset_df)
        mapping = string(r[:MAPPING])
        (ismissing(r[:MAPPING]) || mapping == "NA") && continue
        fg = get(LABELSET_MAPPING_TO_FG, mapping, nothing)
        isnothing(fg) && continue
        desc_dot = string(r[:DESCRIPTION])
        desc_space = replace(desc_dot, "." => " ")
        result[desc_dot] = fg
        result[desc_space] = fg
    end
    return result
end

"""
Build a Dict mapping taxon/species name → Vector of (cscape_group, weight) tuples
from the IPM species mapping file. Entries with dual groups ("X and Y") are split
with weight 0.5 each. `name_col` is the column holding the species name.
"""
function build_taxon_fg_map(
    df::DataFrame, name_col::Symbol
)::Dict{String,Vector{Tuple{String,Float64}}}
    result = Dict{String,Vector{Tuple{String,Float64}}}()
    for r in eachrow(df)
        name = strip(string(r[name_col]))
        g = strip(string(r[:Cscape_group]))
        pairs = if contains(g, " and ")
            [(strip(s), 0.5) for s in split(g, " and ")]
        else
            [(g, 1.0)]
        end
        existing = get(result, name, Tuple{String,Float64}[])
        result[name] = vcat(existing, pairs)
    end
    return result
end

# ─── Load inputs ──────────────────────────────────────────────────────────────

@info "Loading photogrammetry data" file = basename(IPM_FILE)
ipm_photo = CSV.read(IPM_FILE, DataFrame; missingstring=["", "NA", "N/A"])
filter!(row -> row.DATASET == "photogrammetry", ipm_photo)
@info "Photogrammetry rows retained" n = nrow(ipm_photo)

@info "Loading juvenile quadrat data" file = basename(JUV_FILE)
juv_raw = CSV.read(JUV_FILE, DataFrame; missingstring=["", "NA", "N/A"], comment="#")
# Normalise boolean SURVIVAL_USE/GROWTH_USE (true/false) → "yes"/"no" strings
for col in [:SURVIVAL_USE, :GROWTH_USE]
    col ∉ propertynames(juv_raw) && continue
    col_data = juv_raw[!, col]
    eltype(col_data) <: Union{Bool,Missing} || continue
    juv_raw[!, col] = Union{String,Missing}[
        ismissing(v) ? missing : (v ? "yes" : "no") for v in col_data
    ]
end
@info "Juvenile quadrat rows loaded" n = nrow(juv_raw)

ipm = vcat(ipm_photo, juv_raw; cols=:union)

@info "Loading benthic cover data" file = basename(BENTHIC_FILE)
benthic_raw = CSV.read(BENTHIC_FILE, DataFrame; missingstring=["", "NA", "N/A"])

@info "Loading oceanographic stats" file = basename(OCN_STATS_FILE)
ocn = DataFrame(Parquet2.readfile(OCN_STATS_FILE))

@info "Loading species mapping" file = basename(SPECIES_FILE)
species_map = CSV.read(SPECIES_FILE, DataFrame; missingstring=["", "NA"])

@info "Loading benthic labelset" file = basename(LABELSET_FILE)
labelset = CSV.read(LABELSET_FILE, DataFrame; missingstring=[""])

@info "Loading site lookup" file = basename(LOOKUP_FILE)
lookup_file = CSV.read(LOOKUP_FILE, DataFrame; missingstring=["", "NA"])

# ─── Build functional group lookup maps ──────────────────────────────────────

# Benthic photoquadrat: DESCRIPTION → internal functional group (one-to-one, from labelset)
benthic_desc_fg = build_benthic_fg_map(labelset)
@info "Labelset: $(length(benthic_desc_fg) ÷ 2) species mapped to functional groups"

# IPM individual data: taxon → cscape_group (from ecorrap_to_cscape_species.csv)
# photogrammetry dataset rows: match on Updated.name
photo_fg = build_taxon_fg_map(
    species_map[lowercase.(species_map.Dataset) .== "photogrammetry", :],
    :Code
    # Symbol("Updated.name")
)
# juv_quadrat dataset rows: match on Code (genus-level)
juv_fg = build_taxon_fg_map(
    species_map[lowercase.(species_map.Dataset) .== "juv_quadrat", :],
    :Code
)

# ─── Aggregate benthic cover to functional groups ─────────────────────────────

# Identify species cover columns (exclude metadata and _STE standard-error columns)
meta_cols_benthic = Set([:survey_title, :year, :site, :depth_m, :site_reef_name, :transect])
cover_cols = [
    n for n in propertynames(benthic_raw)
    if n ∉ meta_cols_benthic && !endswith(string(n), "_STE")
]

@info "Benthic cover: $(length(cover_cols)) species columns"

# Initialise accumulator DataFrame with one row per (transect × survey_title)
n = nrow(benthic_raw)
parsed_titles = [parse_survey_title(string(t)) for t in benthic_raw.survey_title]
habitat_areas = [p[1] for p in parsed_titles]
depth_cats = [p[2] for p in parsed_titles]

benthic_fg = DataFrame(;
    site_code=string.(benthic_raw.site),
    year=Int.(benthic_raw.year),
    habitat_area=habitat_areas,
    depth_cat=depth_cats
)

for fg in FUNCTIONAL_GROUPS
    benthic_fg[!, Symbol(fg * "_cover")] = zeros(Float64, n)
end

# Distribute each species column to its functional group (one-to-one via labelset)
n_unmapped = 0
for col in cover_cols
    col_str = string(col)
    fg = get(benthic_desc_fg, col_str, nothing)
    if isnothing(fg)
        global n_unmapped += 1
        continue
    end
    out_col = Symbol(fg * "_cover")
    out_col ∈ propertynames(benthic_fg) || continue
    benthic_fg[!, out_col] .+= coalesce.(benthic_raw[!, col], 0.0)
end
n_unmapped > 0 && @warn "Benthic cover columns not found in labelset" n = n_unmapped

# Sum all functional groups to total coral cover per row
benthic_fg.total_coral_cover = [
    sum(benthic_fg[i, Symbol(fg * "_cover")] for fg in FUNCTIONAL_GROUPS)
    for i in 1:nrow(benthic_fg)
]

# Average across transects within (site_code, habitat_area, depth_cat, year)
fg_cover_cols = [Symbol(fg * "_cover") for fg in FUNCTIONAL_GROUPS]
benthic_site = combine(
    groupby(benthic_fg, [:site_code, :habitat_area, :depth_cat, :year]),
    [fg_cover_cols; :total_coral_cover] .=> mean .=> [fg_cover_cols; :total_coral_cover]
)
rename!(benthic_site, :year => :survey_year)

@info "Benthic cover aggregated" rows = nrow(benthic_site) sites = length(
    unique(benthic_site.site_code)
)

# ─── Prepare IPM data ─────────────────────────────────────────────────────────

ipm_work = copy(ipm)

# Standardise column names to lowercase symbols
rename!(ipm_work, [n => Symbol(lowercase(string(n))) for n in propertynames(ipm_work)])

# survey_year from TRANSITION
ipm_work.survey_year = Union{Int,Missing}[
    transition_to_survey_year(string(t)) for t in ipm_work.transition
]

# habitat_area from SITE (Back1 → BA1, etc.)
ipm_work.habitat_area = Union{String,Missing}[
    get(SITE_TO_HABITAT, string(s), missing) for s in ipm_work.site
]

# depth_cat: column is already "S"/"D" in the DEPTH column → rename for join clarity
# (DEPTH has period in name when lowercased; check actual name after rename)
depth_col = :depth
if depth_col ∉ propertynames(ipm_work)
    # Fallback: find any column matching depth
    depth_col = first(filter(c -> occursin("depth", string(c)), propertynames(ipm_work)))
end
rename!(ipm_work, depth_col => :depth_cat)

# Standardise depth_cat to match benthic data (S -> S, D -> D)
# Based on sampling: IPM uses "S" (Shallow) and "D" (Deep).
# Benthic survey_title parsing also yields "S" or "D".
# We ensure they are stripped and uppercase to avoid join failures due to whitespace or case.
ipm_work.depth_cat = [
    ismissing(d) ? missing : uppercase(strip(string(d))) for d in ipm_work.depth_cat
]

# site_code from REEF via lookup (reef names e.g. "moore" → "ONMO")
reef_to_code = Dict(
    zip(
        lowercase.(string.(lookup_file.reef)),
        coalesce.(lookup_file.site_code, "")
    )
)
ipm_work.site_code = Union{String,Missing}[
    let code = get(reef_to_code, lowercase(string(r)), "")
        isempty(code) ? missing : code
    end
    for r in ipm_work.reef
]

# ocn_site_code for oceanographic join (handles Southern GBR naming mismatch)
site_to_ocn = Dict(
    zip(
        coalesce.(lookup_file.site_code, ""),
        coalesce.(lookup_file.ocn_site_code, "")
    )
)
ipm_work.ocn_site_code = Union{String,Missing}[
    let sc = coalesce(ipm_work.site_code[i], "")
        code = get(site_to_ocn, sc, "")
        isempty(code) ? missing : code
    end
    for i in 1:nrow(ipm_work)
]

# Diameter (cm) from tissue area (cm²)
ipm_work.diam = Union{Float64,Missing}[
    ismissing(a) ? missing : Kora.area_to_diam(Float64(a)) for
    a in ipm_work.area_t1_sqcm
]
ipm_work.diamnext = Union{Float64,Missing}[
    ismissing(a) ? missing : Kora.area_to_diam(Float64(a)) for
    a in ipm_work.area_t2_sqcm
]

# Cscape_group from TAXON — choose map based on DATASET
function lookup_cscape(
    dataset_str::AbstractString, taxon_str::AbstractString
)::Union{String,Missing}
    ds = lowercase(strip(dataset_str))
    tx = strip(taxon_str)
    map = ds == "photogrammetry" ? photo_fg : juv_fg
    matches = get(map, tx, Tuple{String,Float64}[])
    isempty(matches) && return missing
    # For dual-group taxa return the first group; the group split is handled in benthic cover
    return matches[1][1]
end

ipm_work.cscape_group = Union{String,Missing}[
    ismissing(ipm_work.taxon[i]) ? missing :
    lookup_cscape(string(ipm_work.dataset[i]), string(ipm_work.taxon[i]))
    for i in 1:nrow(ipm_work)
]

# Extract date_t1 / date_t2 from year-specific date columns (date_2021, date_2022, date_2023)
function get_date(row::DataFrameRow, year::Int)::Union{Date,Missing}
    col = Symbol("date_$(year)")
    col ∉ propertynames(row) && return missing
    v = row[col]
    ismissing(v) && return missing
    v isa Date && return v
    v isa DateTime && return Date(v)
    # Try parsing from string (expected format: DD/MM/YYYY)
    d = tryparse(Date, string(v), dateformat"dd/mm/yyyy")
    isnothing(d) && return missing
    return d
end

date_t1 = Vector{Union{Date,Missing}}(undef, nrow(ipm_work))
date_t2 = Vector{Union{Date,Missing}}(undef, nrow(ipm_work))

for (i, row) in enumerate(eachrow(ipm_work))
    t = string(row.transition)
    parts = split(t, "_")
    if length(parts) == 2
        y1 = tryparse(Int, parts[1])
        y2 = tryparse(Int, parts[2])
        date_t1[i] = isnothing(y1) ? missing : get_date(row, y1)
        date_t2[i] = isnothing(y2) ? missing : get_date(row, y2)
    else
        date_t1[i] = missing
        date_t2[i] = missing
    end
end

ipm_work.date_t1 = date_t1
ipm_work.date_t2 = date_t2

# SURVIVAL_USE: absent → default "yes" (missing colony = genuine mortality, not survey gap)
# GROWTH_USE:   absent → default missing (unknown whether bleaching-affected shrinkage
#               confounds the measurement; downstream analysis should decide how to handle)
if :survival_use ∉ propertynames(ipm_work)
    @warn "survival_use column absent — defaulting to \"yes\" for all rows"
    ipm_work.survival_use = fill("yes", nrow(ipm_work))
end
if :growth_use ∉ propertynames(ipm_work)
    @warn "growth_use column absent — defaulting to missing (bleaching-affected shrinkage unknown)"
    ipm_work.growth_use = Vector{Union{String,Missing}}(missing, nrow(ipm_work))
end

@info "IPM data prepared" rows = nrow(ipm_work) sites = length(
    unique(skipmissing(ipm_work.site_code))
)

# ─── Save aggregated benthic cover ──────────────────────────────────────────

@info "Saving aggregated benthic cover..."
Parquet2.writefile(joinpath(DATA_DIR, "benthic_site.parquet"), benthic_site)

# ─── Join oceanographic stats ─────────────────────────────────────────────────

@info "Joining oceanographic stats..."
ocn_join = select(ocn, Not(:window_type))
rename!(ocn_join, :site_code => :ocn_site_code, :period_year => :survey_year)

expanded = leftjoin(
    ipm_work, ocn_join;
    on=[:ocn_site_code, :depth_cat, :survey_year],
    makeunique=false
)

n_no_ocn = count(row -> ismissing(row.temp_mean_mean), eachrow(expanded))
if n_no_ocn > 0
    @warn "Rows without matching oceanographic data" n = n_no_ocn total = nrow(ipm_work)
    unmatched_ocn = unique(
        expanded[
            ismissing.(expanded.temp_mean_mean),
            [:site_code, :depth_cat, :survey_year]
        ]
    )
    @warn "Unmatched (ocn_site, depth, year) combinations" unmatched_ocn
end

# ─── Select and order output columns ─────────────────────────────────────────

output_cols = [
    # ── Identifiers ───────────────────────────────────────────────────────────
    :dataset, :cluster, :reef, :site_code, :habitat_area, :depth_cat,
    :plot, :colony_id, :taxon, :cscape_group,
    # ── Individual tracking ────────────────────────────────────────────────────
    :transition, :survey_year,
    :diam, :diamnext,
    :area_t1_sqcm, :area_t2_sqcm,
    :survival, :survival_use, :growth_use,
    Symbol("days_t1.t2"),
    :date_t1, :date_t2,
    :bleaching_scores,
    # ── Juvenile-specific ─────────────────────────────────────────────────────
    :quadrat_number, :water_clarity,
    :coral_cover_2021, :coral_cover_2022, :coral_cover_2023,
    # ── Community benthic cover ────────────────────────────────────────────────
    :acro_table_cover, :acro_corym_cover, :corym_non_acro_cover,
    :small_massive_cover, :large_massive_cover, :total_coral_cover,
    # ── Oceanographic ─────────────────────────────────────────────────────────
    :temp_mean_mean, :temp_mean_median,
    :temp_max_mean, :temp_max_median, :n_days_temp,
    :psal_mean_mean, :psal_mean_median, :n_days_psal,
    :cspd_mean_mean, :cspd_mean_median, :n_days_cspd,
    :wave_hs_mean, :wave_hs_median, :n_days_waves,
    :par_dli_mean, :par_dli_median, :n_days_par,
    :depth_min_mean, :depth_min_median,
    :depth_max_mean, :depth_max_median,
    :depth_range_mean, :depth_range_median, :n_days_depth
]

present = Set(propertynames(expanded))
missing_expected = [c for c in output_cols if c ∉ present]
isempty(missing_expected) ||
    @warn "Expected output columns not found" cols = missing_expected

final_cols = [c for c in output_cols if c ∈ present]
output = expanded[:, final_cols]
sort!(output, :survey_year)

@info "Writing expanded dataset" rows = nrow(output) cols = ncol(output) file = OUTPUT_FILE
mkpath(dirname(OUTPUT_FILE))
Parquet2.writefile(OUTPUT_FILE, output)

@info "Done. Output: $OUTPUT_FILE"
