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

# Write to file
write(REPORT_FILE, String(take!(io)))
@info "Report written" file = REPORT_FILE
