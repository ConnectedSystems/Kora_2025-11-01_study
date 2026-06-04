"""
Generate a markdown coverage report for ecorrap_expanded.parquet.

Reports, for each column:
  - total missing/NaN rows
  - breakdown by survey_year

Output: data/coverage_report.md
"""

using CSV
using DataFrames
using Parquet2

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const INPUT_FILE = joinpath(DATA_DIR, "ecorrap_expanded.parquet")
const REPORT_FILE = joinpath(DATA_DIR, "coverage_report.md")



# ─── Helpers ──────────────────────────────────────────────────────────────────

is_absent(v) = ismissing(v) || (v isa AbstractFloat && isnan(v))

function coverage_by_year(df::DataFrame, col::Symbol)
    years = sort(unique(skipmissing(df.survey_year)))
    rows = NamedTuple{(:survey_year, :n_total, :n_present, :n_missing, :n_nan, :pct_present),
                      Tuple{Int,Int,Int,Int,Int,Float64}}[]
    for yr in years
        sub = df[df.survey_year .== yr, :]
        n = nrow(sub)
        vals = sub[!, col]
        n_nan  = count(v -> !ismissing(v) && v isa AbstractFloat && isnan(v), vals)
        n_miss = count(ismissing, vals)
        n_pres = n - n_nan - n_miss
        push!(rows, (survey_year=yr, n_total=n, n_present=n_pres,
                     n_missing=n_miss, n_nan=n_nan,
                     pct_present=round(100 * n_pres / n; digits=1)))
    end
    return rows
end

# ─── Load ─────────────────────────────────────────────────────────────────────

@info "Loading expanded dataset" file = INPUT_FILE
df = DataFrame(Parquet2.readfile(INPUT_FILE))
@info "Loaded" rows = nrow(df) cols = ncol(df)

years_all = sort(unique(skipmissing(df.survey_year)))
datasets  = sort(unique(skipmissing(df.dataset)))

# ─── Column groups ────────────────────────────────────────────────────────────

id_cols = [:dataset, :cluster, :reef, :site_code, :habitat_area, :depth_cat,
           :plot, :colony_id, :taxon, :cscape_group]

tracking_cols = [:transition, :survey_year, :diam, :diamnext,
                 :area_t1_sqcm, :area_t2_sqcm, :survival, :survival_use, :growth_use,
                 Symbol("days_t1.t2"), :date_t1, :date_t2, :bleaching_scores]

juv_cols = [:quadrat_number, :water_clarity,
            :coral_cover_2021, :coral_cover_2022, :coral_cover_2023]

ocn_cols = [:temp_mean_mean, :temp_mean_median, :temp_max_mean, :temp_max_median, :n_days_temp,
            :psal_mean_mean, :psal_mean_median, :n_days_psal,
            :cspd_mean_mean, :cspd_mean_median, :n_days_cspd,
            :wave_hs_mean, :wave_hs_median, :n_days_waves,
            :par_dli_mean, :par_dli_median, :n_days_par]

# ─── Build report ─────────────────────────────────────────────────────────────

io = IOBuffer()

println(io, "# EcoRRAP Expanded Dataset — Data Coverage Report")
println(io, "")
println(io, "Generated from: `$(basename(INPUT_FILE))`  ")
println(io, "Total rows: **$(nrow(df))**  ")
println(io, "Survey years: $(join(years_all, ", "))  ")
println(io, "Datasets: $(join(datasets, ", "))")
println(io, "")

# Row counts per dataset × year
println(io, "## Row counts by dataset and survey year")
println(io, "")
println(io, "| dataset | " * join(string.(years_all), " | ") * " | **total** |")
println(io, "|---------|" * repeat("--------|", length(years_all)) * "----------|")
for ds in datasets
    sub = df[df.dataset .=== ds, :]
    counts = [count(==(yr), sub.survey_year) for yr in years_all]
    println(io, "| $ds | " * join(string.(counts), " | ") * " | **$(nrow(sub))** |")
end
println(io, "")

# Per-column coverage table helper
function section(io, title, cols)
    println(io, "## $title")
    println(io, "")
    println(io, "Rows with usable data (non-missing, non-NaN) per survey year.  ")
    println(io, "Columns absent from the dataset are skipped.")
    println(io, "")

    present_cols = filter(c -> c ∈ propertynames(df), cols)

    for col in present_cols
        rows = coverage_by_year(df, col)
        isempty(rows) && continue

        println(io, "### `$col`")
        println(io, "")
        println(io, "| survey_year | n_total | present | missing | NaN | % present |")
        println(io, "|-------------|---------|---------|---------|-----|-----------|")
        for r in rows
            println(io, "| $(r.survey_year) | $(r.n_total) | $(r.n_present) | $(r.n_missing) | $(r.n_nan) | $(r.pct_present)% |")
        end

        total_present = sum(r.n_present for r in rows)
        total_rows    = sum(r.n_total for r in rows)
        pct_total     = round(100 * total_present / total_rows; digits=1)
        println(io, "| **total** | **$total_rows** | **$total_present** | | | **$pct_total%** |")
        println(io, "")
    end
end

section(io, "Identifier and taxonomic columns", id_cols)
section(io, "Individual tracking columns", tracking_cols)
section(io, "Juvenile-specific columns", juv_cols)
section(io, "Oceanographic columns", ocn_cols)

# NaN-specific note
println(io, "---")
println(io, "")
println(io, "## Notes on NaN vs missing")
println(io, "")
println(io, """
Empty cells in the output represent rows where no data existed for the join key
(site × depth × year). This is expected — not all sites have all data types.

`NaN` values in oceanographic columns arise when a sensor recorded data for a period
but all individual measurements were flagged as bad quality (so `skipmissing()` returns
an empty collection). These are treated identically to missing by downstream analyses
and were replaced with `missing` in the aggregation pipeline (fix applied in
`00a_prep_oceanographic_stats.jl`).

Affected combinations (pre-fix):

| variable | site_code | depth_cat | survey_year |
|----------|-----------|-----------|-------------|
| psal     | OSHE      | D         | 2021        |
| par      | OSKE      | S         | 2022        |
| par      | TSMA      | D         | 2022        |
""")

# ─── Subsets for study variable summary (expanded parquet only) ──────────────

on_df = df[coalesce.(df.cluster .== "offshore_north", false), :]
ts_df = df[coalesce.(df.cluster .== "torres_strait",  false), :]

on_g = on_df[coalesce.(on_df.growth_use   .== "yes", false), :]
on_s = on_df[coalesce.(on_df.survival_use .== "yes", false), :]
ts_g = ts_df[coalesce.(ts_df.growth_use   .== "yes", false), :]
ts_s = ts_df[coalesce.(ts_df.survival_use .== "yes", false), :]

# ─── Study variable summary table ────────────────────────────────────────────

# Helpers
_vals(d, c)  = collect(skipmissing(d[!, c]))
_uvals(d, c) = sort(unique(skipmissing(d[!, c])))
_nu(d, c)    = length(_uvals(d, c))

function _fmt_cr(vals::Vector)   # count / [range]
    isempty(vals) && return "N/A"
    lo = round(minimum(vals); digits=2)
    hi = round(maximum(vals); digits=2)
    return "$(length(vals)) / [$(lo) – $(hi)]"
end

function _site_stat(df, site_col::Symbol, val_col::Symbol)
    rows = unique(dropmissing(df[:, [site_col, val_col]]))
    vals = collect(skipmissing(rows[!, val_col]))
    return _fmt_cr(vals)
end

println(io, "---")
println(io, "")
println(io, "## Study variable summary")
println(io, "")
println(io, """
Counts and value ranges for key variables in the Offshore North and Torres Strait
regional analysis datasets.  Diameter counts are filtered to `growth_use = "yes"` /
`survival_use = "yes"`.  All values are sourced from `ecorrap_expanded.parquet`.
Continuous depth for Torres Strait uses `depth_min_mean` from the EcoRRAP in-situ
logger deployment at Masig Reef; only Masig has logger coverage in the current
study scope.
""")

println(io, "| Identifier | Description | Offshore North<br>Count / [Range] | Torres Strait<br>Count / [Range] |")
println(io, "|------------|-------------|-----------------------------------|----------------------------------|")

# diam / diam_mort  (diam = area_to_diam(area_t1), diam_mort = area_to_diam(area_t2) for survival)
on_g_diam = _vals(on_g, :diam)
on_s_diam = _vals(on_s, :diamnext)
ts_g_diam = _vals(ts_g, :diam)
ts_s_diam = _vals(ts_s, :diamnext)
on_diam_cell = "Growth: $(_fmt_cr(on_g_diam))<br>Survival: $(_fmt_cr(on_s_diam))"
ts_diam_cell = "Growth: $(_fmt_cr(ts_g_diam))<br>Survival: $(_fmt_cr(ts_s_diam))"
println(io, "| `diam` / `diam_mort` | Estimated coral diameter at observation (cm): `diam` for growth obs., `diam_mort` for mortality obs. | $(on_diam_cell) | $(ts_diam_cell) |")

# plot / quadrat  (unique plot values per cluster)
on_plot_n = _nu(on_g, :plot)
ts_quad_n = _nu(ts_g, :quadrat_number)
println(io, "| `plot` / `quadrat` | Specific monitored plot (ON) or quadrat (TS) | $(on_plot_n) | $(ts_quad_n) |")

# depth — depth_cat (S/D) is available for both; continuous depth_min_mean from logger for TS only
on_depth_str = join(sort(unique(skipmissing(on_g.depth_cat))), ", ")
ts_depth_str = _site_stat(
    unique(dropmissing(ts_df[:, [:site_code, :depth_cat, :depth_min_mean]])),
    :site_code, :depth_min_mean
)
println(io, "| `depth` | Depth category (S/D) for ON; continuous logger mean (m) for TS (Masig only) | $(on_depth_str) | $(ts_depth_str) |")

# temperature — EcoRRAP in-situ logger data for both regions
# ON: temp_max_mean (mean of daily maxima); TS: temp_mean_mean (mean of daily means)
on_temp_str = _site_stat(
    unique(dropmissing(on_df[:, [:site_code, :depth_cat, :temp_max_mean]])),
    :site_code, :temp_max_mean
)
ts_temp_str = _site_stat(
    unique(dropmissing(ts_df[:, [:site_code, :depth_cat, :temp_mean_mean]])),
    :site_code, :temp_mean_mean
)
println(io, "| `temp` | Mean of daily temperature (°C) from EcoRRAP in-situ loggers; daily-max mean for ON, daily mean for TS | $(on_temp_str) | $(ts_temp_str) |")

# habitat
on_hab = _uvals(on_g, :habitat_area)
ts_hab = _uvals(ts_g, :habitat_area)
on_hab_str = "$(length(on_hab)) ($(join(on_hab, "; ")))"
ts_hab_str = "$(length(ts_hab)) ($(join(ts_hab, "; ")))"
println(io, "| `habitat` | Reef habitat zone code (`habitat_area`) | $(on_hab_str) | $(ts_hab_str) |")

# long / lat (site-level; from benthic CSV metadata via 00b)
on_lon_str = _site_stat(unique(dropmissing(on_df[:, [:site_code, :lon]])), :site_code, :lon)
on_lat_str = _site_stat(unique(dropmissing(on_df[:, [:site_code, :lat]])), :site_code, :lat)
ts_lon_str = _site_stat(unique(dropmissing(ts_df[:, [:site_code, :lon]])), :site_code, :lon)
ts_lat_str = _site_stat(unique(dropmissing(ts_df[:, [:site_code, :lat]])), :site_code, :lat)
println(io, "| `lon` | Longitude (decimal degrees) | $(on_lon_str) | $(ts_lon_str) |")
println(io, "| `lat` | Latitude (decimal degrees) | $(on_lat_str) | $(ts_lat_str) |")

# site (reef-level code in expanded parquet)
on_site_n = _nu(on_g, :site_code)
ts_site_n = _nu(ts_g, :site_code)
println(io, "| `site` | Reef-level site code | $(on_site_n) | $(ts_site_n) |")

# reef
on_reef = _uvals(on_g, :reef)
ts_reef = _uvals(ts_g, :reef)
on_reef_str = "$(length(on_reef)) ($(join(titlecase.(on_reef), "; ")))"
ts_reef_str = "$(length(ts_reef)) ($(join(titlecase.(ts_reef), "; ")))"
println(io, "| `reef` | Monitored reef | $(on_reef_str) | $(ts_reef_str) |")

# taxa
on_g_taxa = _nu(on_g, :taxon)
on_s_taxa = _nu(on_s, :taxon)
ts_g_taxa = _nu(ts_g, :taxon)
ts_s_taxa = _nu(ts_s, :taxon)
println(io, "| `taxa` | Individual taxonomic group | Growth: $(on_g_taxa)<br>Survival: $(on_s_taxa) | Growth: $(ts_g_taxa)<br>Survival: $(ts_s_taxa) |")

# functional_group
on_fg_n = _nu(on_g, :cscape_group)
ts_fg_n = _nu(ts_g, :cscape_group)
println(io, "| `functional_group` | Morphological grouping of taxa | $(on_fg_n) | $(ts_fg_n) |")

# colony_id
ts_g_cid = _nu(ts_g, :colony_id)
ts_s_cid = _nu(ts_s, :colony_id)
println(io, "| `colony_id` | Individual coral identifier | N/A | Growth: $(ts_g_cid)<br>Survival: $(ts_s_cid) |")
println(io, "")

# ─── Habitat type codes ───────────────────────────────────────────────────────

habitat_zone_labels = Dict("BA" => "Back reef", "FL" => "Flank", "FR" => "Front reef", "LA" => "Lagoon")

println(io, "---")
println(io, "")
println(io, "## Habitat type codes")
println(io, "")
println(io, """
The `habitat_area` column identifies the reef zone of a survey plot or quadrat.
Codes use a two-letter zone prefix followed by a numeric plot index
(e.g. `BA1` = first back-reef plot, `FR2` = second front-reef plot).
The `habitat` column in model inputs uses the zone label directly
(`back`, `flank`, `front`, `lagoon`).
""")
println(io, "| Prefix | Zone | Present in |")
println(io, "|--------|------|------------|")
println(io, "| `BA` | Back reef | Offshore North, Torres Strait |")
println(io, "| `FL` | Flank | Offshore North, Torres Strait |")
println(io, "| `FR` | Front reef | Offshore North, Torres Strait |")
println(io, "| `LA` | Lagoon | Torres Strait only |")
println(io, "")

all_hab_codes = sort(unique(skipmissing(df.habitat_area)))
println(io, "Full `habitat_area` inventory in the expanded dataset:")
println(io, "")
println(io, "| Code | Zone | Plot index | Datasets present |")
println(io, "|------|------|------------|------------------|")
for code in all_hab_codes
    prefix = code[1:2]
    idx    = length(code) >= 3 ? code[3:end] : ""
    zone   = get(habitat_zone_labels, prefix, "Unknown")
    ds_mask = coalesce.(df.habitat_area .== code, false)
    ds_list = sort(unique(skipmissing(df[ds_mask, :cluster])))
    println(io, "| `$(code)` | $(zone) | $(idx) | $(join(ds_list, ", ")) |")
end
println(io, "")

# Write to file
write(REPORT_FILE, String(take!(io)))
@info "Report written" file = REPORT_FILE
