"""
Generate a paper-ready dataset summary from the fully-filtered pipeline CSVs.

These CSVs are the outputs of the 01b / 01c exploratory SA scripts and reflect the
exact observations passed to every downstream script (02 and 03 series):
  - lin_ext > 0 filter applied (partial-mortality / shrinking colonies excluded)
  - only taxa that map to one of the 5 TARGET_GROUPS functional groups retained

Reads:
  data/offshore_north/overall/offshore_north_growth.csv
  data/offshore_north/overall/offshore_north_survival.csv
  data/torres_strait/overall/torres_strait_growth.csv
  data/torres_strait/overall/torres_strait_survival.csv

Outputs:
  data/logger_report/paper_dataset_summary.md          — full diagnostic summary
  data/logger_report/paper_table1_dataset_overview.md  — Table 1 body, included by paper.qmd
"""

using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "common.jl"))

let

REPORT_FILE = joinpath(OUTPUT_DIR, "logger_report", "paper_dataset_summary.md")
TABLE1_FILE = joinpath(OUTPUT_DIR, "logger_report", "paper_table1_dataset_overview.md")

# ─── Load pipeline output CSVs ────────────────────────────────────────────────

on_g = CSV.read(
    joinpath(OUTPUT_DIR, "offshore_north", "overall", "offshore_north_growth.csv"),
    DataFrame
)
on_s = CSV.read(
    joinpath(OUTPUT_DIR, "offshore_north", "overall", "offshore_north_survival.csv"),
    DataFrame
)
ts_g = CSV.read(
    joinpath(OUTPUT_DIR, "torres_strait", "overall", "torres_strait_growth.csv"),
    DataFrame
)
ts_s = CSV.read(
    joinpath(OUTPUT_DIR, "torres_strait", "overall", "torres_strait_survival.csv"),
    DataFrame
)

@info "Loaded pipeline CSVs" on_g=nrow(on_g) on_s=nrow(on_s) ts_g=nrow(ts_g) ts_s=nrow(ts_s)

# ─── Helpers ──────────────────────────────────────────────────────────────────

day_col = Symbol("days_t1.t2")

"""Interval statistics excluding 0-day entries (same-survey artefacts)."""
function interval_stats(df::DataFrame)
    days = filter(d -> d > 0, collect(skipmissing(df[:, day_col])))
    isempty(days) && return (mean=missing, sd=missing, min=missing, max=missing)
    return (
        mean = round(mean(days); digits=1),
        sd   = round(std(days);  digits=1),
        min  = minimum(days),
        max  = maximum(days),
    )
end

"""Unique (reef × habitat_type) combinations — the paper's 'monitoring sites'."""
n_sites(df) = nrow(unique(df[:, [:reef, :site_code, :habitat_type]]))

"""Unique (reef × habitat_type × depth_cat × plot) combinations — the paper's 'plots'."""
n_plots(df) = nrow(unique(df[:, [:reef, :site_code, :habitat_type, :depth_cat, :plot]]))

hab_zones(df) = sort(unique(skipmissing(df.habitat_type)))

ZONE_LABEL = Dict("back" => "back", "flank" => "flank", "front" => "front", "lagoon" => "lagoon")

function zone_list(df)
    zones = hab_zones(df)
    labels = [get(ZONE_LABEL, z, z) for z in zones]
    return "$(length(zones)) ($(join(labels, ", ")))"
end

reefs_str(df) = join(titlecase.(sort(unique(skipmissing(df.reef)))), ", ")
taxa_n(df)    = length(unique(skipmissing(df.taxa)))
transitions_str(df) = join(sort(unique(skipmissing(df.transition))), ", ")

# ─── Compute per-region stats ─────────────────────────────────────────────────

function region_stats(g::DataFrame, s::DataFrame)
    g_int = interval_stats(g)
    s_int = interval_stats(s)
    return (
        growth_n      = nrow(g),
        survival_n    = nrow(s),
        plots         = n_plots(g),
        sites         = n_sites(g),
        taxa_g        = taxa_n(g),
        taxa_s        = taxa_n(s),
        habitats      = zone_list(g),
        reefs         = reefs_str(g),
        transitions   = transitions_str(g),
        g_mean        = g_int.mean,
        g_sd          = g_int.sd,
        g_min         = g_int.min,
        g_max         = g_int.max,
        s_mean        = s_int.mean,
        s_sd          = s_int.sd,
        s_min         = s_int.min,
        s_max         = s_int.max,
    )
end

on = region_stats(on_g, on_s)
ts = region_stats(ts_g, ts_s)

# ─── Write report ─────────────────────────────────────────────────────────────

io = IOBuffer()

println(io, """
# Paper dataset summary

Numbers derived from the fully-filtered pipeline CSVs
(`offshore_north/overall/` and `torres_strait/overall/`).
Filters applied on top of `ecorrap_expanded.parquet`:
1. `lin_ext > 0` — excludes partial-mortality / shrinking colonies (growth only)
2. Taxa filter — only taxa mapping to one of the 5 TARGET_GROUPS functional groups retained
""")

# ─── Per-region summary tables ────────────────────────────────────────────────

for (label, st) in [("Offshore North", on), ("Torres Strait", ts)]
    println(io, """
    ---

    ## $label

    | Statistic | Value |
    |-----------|-------|
    | Reefs | $(st.reefs) |
    | Transitions (survey years) | $(st.transitions) |
    | Growth observations | $(st.growth_n) |
    | Survival observations | $(st.survival_n) |
    | Monitoring sites (reef × habitat_type) | $(st.sites) |
    | Plots (reef × habitat_type × plot) | $(st.plots) |
    | Taxa — growth | $(st.taxa_g) |
    | Taxa — survival | $(st.taxa_s) |
    | Habitat types | $(st.habitats) |
    | Mean survey interval — growth (days) | $(st.g_mean) (SD: $(st.g_sd); range: $(st.g_min)–$(st.g_max)) |
    | Mean survey interval — survival (days) | $(st.s_mean) (SD: $(st.s_sd); range: $(st.s_min)–$(st.s_max)) |

    """)
end

# ─── Per-transition breakdown ─────────────────────────────────────────────────

println(io, """
---

## Observation counts by reef and transition

""")

function print_transition_table(io::IOBuffer, label::String, kind::String, df::DataFrame)
    println(io, "### $label — $kind\n")
    counts = sort(combine(groupby(df, [:reef, :transition]), nrow => :n), [:reef, :transition])
    println(io, "| reef | transition | n |\n|------|------------|---|")
    for r in eachrow(counts)
        println(io, "| $(r.reef) | $(r.transition) | $(r.n) |")
    end
    println(io, "| **total** | | **$(nrow(df))** |\n")
end

for (label, g, s) in [("Offshore North", on_g, on_s), ("Torres Strait", ts_g, ts_s)]
    print_transition_table(io, label, "growth", g)
    print_transition_table(io, label, "survival", s)
end

# ─── Variable summary table (matches coverage_report.md style) ───────────────

println(io, """
---

## Study variable summary (pipeline-filtered)

Sourced from the pipeline CSVs (post lin_ext and taxa filters).
Diameter ranges from `diam` (growth obs.) and `diam_mort` (survival obs.).
Depth: both regions use plot-level survey depth (`depth_cont`) from the IPM photogrammetry
file, now available for all reefs including Masig. Torres Strait also has continuous logger
mean depth (`depth_min_mean`) from EcoRRAP in-situ loggers (available at Masig Reef only) for
comparison.
Temperature: mean of daily maxima for ON (`temp_max_mean`), mean of daily means for TS
(`temp_mean_mean`). Colony ID is only tracked in the Torres Strait photogrammetry dataset.
""")

_nu(df, c)   = length(unique(skipmissing(df[!, c])))
_vals(df, c) = collect(skipmissing(df[!, c]))

function _fmt_cr(vals)
    isempty(vals) && return "N/A"
    lo = round(minimum(vals); digits=2)
    hi = round(maximum(vals); digits=2)
    return "$(length(vals)) / [$(lo) – $(hi)]"
end

function _site_stat(df, site_col, val_col)
    rows = unique(dropmissing(df[:, [site_col, val_col]]))
    return _fmt_cr(collect(skipmissing(rows[!, val_col])))
end

function _hab_str(df)
    codes = sort(unique(skipmissing(df.habitat_type)))
    return "$(length(codes)) ($(join(codes, "; ")))"
end

function _reef_str(df)
    reefs = sort(unique(skipmissing(df.reef)))
    return "$(length(reefs)) ($(join(titlecase.(reefs), "; ")))"
end

# diam / diam_mort
on_g_diam = _fmt_cr(_vals(on_g, :diam))
on_s_diam = _fmt_cr(_vals(on_s, :diam_mort))
ts_g_diam = _fmt_cr(_vals(ts_g, :diam))
ts_s_diam = _fmt_cr(_vals(ts_s, :diam_mort))

# plots — composite physical plot identity (reef × site_code × habitat_type × depth_cat ×
# plot/quadrat label), matching the narrative's plot counts and the summary table above
# (on.plots / ts.plots from region_stats). TS quadrats nest cleanly within this composite
# plot definition, so both regions are now reported on the same "plot" basis.
on_plot_n = on.plots
ts_plot_n = ts.plots

# depth — ON: plot-level survey depth (`depth_cont`); TS: continuous logger mean depth
# (`depth_min_mean`, available at Masig Reef only). Categorical `depth_cat` is not used here.
function _depth_stat(g, s, depth_col::Symbol)
    rows = unique(vcat(
        dropmissing(g[:, [:site_code, :habitat_type, depth_col]]),
        dropmissing(s[:, [:site_code, :habitat_type, depth_col]])
    ))
    return _fmt_cr(rows[!, depth_col])
end

on_depth_str = _depth_stat(on_g, on_s, :depth_cont)
ts_depth_str = _depth_stat(ts_g, ts_s, :depth_min_mean)

# temperature (site-level unique values / range)
on_temp_str = _site_stat(
    unique(dropmissing(on_g[:, [:site_code, :depth_cat, :temp_max_mean]])),
    :site_code, :temp_max_mean
)
ts_temp_str = _site_stat(
    unique(dropmissing(ts_g[:, [:site_code, :depth_cat, :temp_mean_mean]])),
    :site_code, :temp_mean_mean
)

# habitat_type
on_hab_str = _hab_str(on_g)
ts_hab_str = _hab_str(ts_g)

# lon / lat (plot-level; renamed from bare :lon/:lat after data update)
on_lon_str = _site_stat(unique(dropmissing(on_g[:, [:site_code, :plot_lon]])), :site_code, :plot_lon)
on_lat_str = _site_stat(unique(dropmissing(on_g[:, [:site_code, :plot_lat]])), :site_code, :plot_lat)
ts_lon_str = _site_stat(unique(dropmissing(ts_g[:, [:site_code, :plot_lon]])), :site_code, :plot_lon)
ts_lat_str = _site_stat(unique(dropmissing(ts_g[:, [:site_code, :plot_lat]])), :site_code, :plot_lat)

# site / reef
on_site_n = _nu(on_g, :site_code)
ts_site_n = _nu(ts_g, :site_code)
on_reef_str = _reef_str(on_g)
ts_reef_str = _reef_str(ts_g)

# taxa (from Cscape_group-filtered data, so uses :taxa column)
on_taxa_g_n = _nu(on_g, :taxa)
on_taxa_s_n = _nu(on_s, :taxa)
ts_taxa_g_n = _nu(ts_g, :taxa)
ts_taxa_s_n = _nu(ts_s, :taxa)

# functional group (Cscape_group)
on_fg_n = _nu(on_g, :Cscape_group)
ts_fg_n = _nu(ts_g, :Cscape_group)

# colony_id (TS photogrammetry only; ON plot-based so N/A)
on_cid = "N/A"
ts_cid_g = _nu(ts_g, :colony_id)
ts_cid_s = _nu(ts_s, :colony_id)

println(io, """
| Identifier | Description | Offshore North<br>Count / [Range] | Torres Strait<br>Count / [Range] |
|------------|-------------|-----------------------------------|----------------------------------|
| `diam` / `diam_mort` | Estimated coral diameter at observation (cm): `diam` for growth obs., `diam_mort` for mortality obs. | Growth: $(on_g_diam)<br>Survival: $(on_s_diam) | Growth: $(ts_g_diam)<br>Survival: $(ts_s_diam) |
| `plot` | Physical monitoring plot (reef × site × habitat × depth × plot/quadrat label) | $(on_plot_n) | $(ts_plot_n) |
| `depth` | Continuous depth (m): ON uses plot-level survey depth (`depth_cont`) from IPM photogrammetry file; TS uses logger mean depth (`depth_min_mean`) from EcoRRAP in-situ loggers (Masig Reef only) | $(on_depth_str) | $(ts_depth_str) |
| `temp` | Mean of daily temperature (°C) from EcoRRAP in-situ loggers; daily-max mean for ON, daily mean for TS | $(on_temp_str) | $(ts_temp_str) |
| `habitat_type` | Reef exposure category (`habitat_type`, e.g. "back"/"front") | $(on_hab_str) | $(ts_hab_str) |
| `lon` | Longitude (decimal degrees) | $(on_lon_str) | $(ts_lon_str) |
| `lat` | Latitude (decimal degrees) | $(on_lat_str) | $(ts_lat_str) |
| `site` | Reef-level site code | $(on_site_n) | $(ts_site_n) |
| `reef` | Monitored reef | $(on_reef_str) | $(ts_reef_str) |
| `taxa` | Individual taxonomic group | Growth: $(on_taxa_g_n)<br>Survival: $(on_taxa_s_n) | Growth: $(ts_taxa_g_n)<br>Survival: $(ts_taxa_s_n) |
| `functional_group` | Morphological grouping of taxa | $(on_fg_n) | $(ts_fg_n) |
| `colony_id` | Individual coral identifier | $(on_cid) | Growth: $(ts_cid_g)<br>Survival: $(ts_cid_s) |
""")

# ─── Table 1 (paper.qmd include) ──────────────────────────────────────────────
#
# Body of "Table 1. Overview of factors assessed in the sensitivity analysis".
# The 12 rows must stay in one-to-one correspondence with the factors in
# `growth_include_cols` / `surv_include_cols` in common.jl — if a factor is added to
# or removed from the SA, this list changes with it.
#
# Cell convention (mirrors the table caption):
#   - `diam`/`diam_mort`: observation count and range over all observations.
#   - `taxa`: unique count, reported separately for growth and survival.
#   - everything else: unique-value count over the growth dataset, with the range
#     across those unique values where the factor is continuous.
#   - a factor with no non-missing values in a region is reported as N/A.
#   - `**` marks a factor that is site-invariant within the region (one unique
#     value) and therefore contributes near-zero PAWN sensitivity.

_present(df, c) = string(c) in names(df) && !isempty(collect(skipmissing(df[!, c])))
_uvals(df, c)   = unique(collect(skipmissing(df[!, c])))

"""Unique-value count, plus [min – max] across those values for continuous factors."""
function t1_unique(df::DataFrame, c::Symbol; with_range::Bool=true)
    _present(df, c) || return ("N/A", false)
    v = _uvals(df, c)
    invariant = length(v) == 1
    with_range || return ("$(length(v))", invariant)
    lo = round(minimum(v); digits=2)
    hi = round(maximum(v); digits=2)
    return ("$(length(v)) / [$(lo) – $(hi)]", invariant)
end

"""Observation count plus value range, for the diameter columns."""
function t1_obs(df::DataFrame, c::Symbol)
    _present(df, c) || return "N/A"
    v = collect(skipmissing(df[!, c]))
    return "$(length(v)) / [$(round(minimum(v); digits=2)) – $(round(maximum(v); digits=2))]"
end

"""Zone count with human-readable labels, e.g. `4 (back; flank; front; lagoon)`."""
function t1_habitat_type(df::DataFrame)
    _present(df, :habitat_type) || return "N/A"
    z = sort(_uvals(df, :habitat_type))
    return "$(length(z)) ($(join(z, "; ")))"
end

# Rows that report a single value shared by both growth and survival are computed
# from the growth dataset; the diameter and taxa rows report both explicitly.
struct T1Row
    id::String
    description::String
    source::String
    on::String
    ts::String
end

"""Build one row, appending the site-invariance marker to the description if earned."""
function t1_row(id, desc, source, (on_cell, on_inv), (ts_cell, ts_inv))
    marker = (on_inv || ts_inv) ? " **" : ""
    return T1Row(id, desc * marker, source, on_cell, ts_cell)
end

table1_rows = T1Row[
    T1Row(
        "`diam` / `diam_mort`",
        "Estimated coral diameter at observation (cm): `diam` for growth obs., " *
        "`diam_mort` for mortality obs.",
        "EcoRRAP photogrammetry colony tracking",
        "Growth: $(t1_obs(on_g, :diam))<br>Survival: $(t1_obs(on_s, :diam_mort))",
        "Growth: $(t1_obs(ts_g, :diam))<br>Survival: $(t1_obs(ts_s, :diam_mort))",
    ),
    T1Row(
        "`taxa`",
        "Individual taxonomic group",
        "EcoRRAP taxonomic colony labels",
        "Growth: $(taxa_n(on_g))<br>Survival: $(taxa_n(on_s))",
        "Growth: $(taxa_n(ts_g))<br>Survival: $(taxa_n(ts_s))",
    ),
    t1_row(
        "`functional_group`", "Morphological grouping of taxa",
        "Kora group mapping from taxa",
        t1_unique(on_g, :functional_group; with_range=false),
        t1_unique(ts_g, :functional_group; with_range=false),
    ),
    T1Row(
        "`habitat_type`", "Reef exposure category",
        "EcoRRAP site/habitat annotations",
        t1_habitat_type(on_g), t1_habitat_type(ts_g),
    ),
    t1_row(
        "`quadrat_number`", "Monitored plot (ON) or quadrat (TS) identifier",
        "EcoRRAP monitoring design metadata",
        t1_unique(on_g, :quadrat_number; with_range=false),
        t1_unique(ts_g, :quadrat_number; with_range=false),
    ),
    t1_row(
        "`depth_cont`", "Plot-level survey depth (m) from IPM photogrammetry",
        "EcoRRAP IPM photogrammetry metadata",
        t1_unique(on_g, :depth_cont), t1_unique(ts_g, :depth_cont),
    ),
    t1_row(
        "`ereefs_temp_mean`",
        "Mean of daily mean sea water temperature (°C) from eReefs GBR1 Hydro v2 " *
        "(~1 km resolution) at −9 m depth over EcoRRAP survey period",
        "eReefs hydrodynamic model",
        t1_unique(on_g, :ereefs_temp_mean), t1_unique(ts_g, :ereefs_temp_mean),
    ),
    t1_row(
        "`ereefs_temp_max`",
        "Mean of daily maximum sea water temperature (°C) from eReefs GBR1 Hydro v2 " *
        "(~1 km resolution) at −9 m depth over EcoRRAP survey period",
        "eReefs hydrodynamic model",
        t1_unique(on_g, :ereefs_temp_max), t1_unique(ts_g, :ereefs_temp_max),
    ),
    t1_row(
        "`wave_ubed90`", "Wave-induced bottom current speed (m/s; 90th percentile)",
        "eReefs hydrodynamic model",
        t1_unique(on_g, :wave_ubed90), t1_unique(ts_g, :wave_ubed90),
    ),
    t1_row(
        "`depth_bathy_m`", "Modelled bathymetric depth (m)",
        "eReefs hydrodynamic model",
        t1_unique(on_g, :depth_bathy_m), t1_unique(ts_g, :depth_bathy_m),
    ),
    t1_row(
        "`plot_lon`", "Longitude (decimal degrees)",
        "EcoRRAP site geospatial metadata",
        t1_unique(on_g, :plot_lon), t1_unique(ts_g, :plot_lon),
    ),
    t1_row(
        "`plot_lat`", "Latitude (decimal degrees)",
        "EcoRRAP site geospatial metadata",
        t1_unique(on_g, :plot_lat), t1_unique(ts_g, :plot_lat),
    ),
]

t1 = IOBuffer()
println(t1, "| Identifier | Description | Data source | Offshore North<br>Count / [Range] | Torres Strait<br>Count / [Range] |")
println(t1, "|------------|-------------|-------------|-----------------------------------|----------------------------------|")
for r in table1_rows
    println(t1, "| $(r.id) | $(r.description) | $(r.source) | $(r.on) | $(r.ts) |")
end

# ─── Write to file ────────────────────────────────────────────────────────────

mkpath(dirname(REPORT_FILE))
write(REPORT_FILE, String(take!(io)))
@info "Report written" file = REPORT_FILE

write(TABLE1_FILE, String(take!(t1)))
@info "Table 1 written" file = TABLE1_FILE n_factors = length(table1_rows)

end # let
