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

Output: data/logger_report/paper_dataset_summary.md
"""

using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "common.jl"))

let

REPORT_FILE = joinpath(OUTPUT_DIR, "logger_report", "paper_dataset_summary.md")

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

"""Unique (reef × habitat_area) combinations — the paper's 'monitoring sites'."""
n_sites(df) = nrow(unique(df[:, [:reef, :site_code, :habitat_area]]))

"""Unique (reef × habitat_area × depth_cat × plot) combinations — the paper's 'plots'."""
n_plots(df) = nrow(unique(df[:, [:reef, :site_code, :habitat_area, :depth_cat, :plot]]))

hab_zones(df) = sort(unique(map(h -> h[1:2], unique(skipmissing(df.habitat_area)))))

ZONE_LABEL = Dict("BA" => "back", "FL" => "flank", "FR" => "front", "LA" => "lagoon")

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
    | Monitoring sites (reef × habitat) | $(st.sites) |
    | Plots (reef × habitat × plot) | $(st.plots) |
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
Depth: Offshore North uses plot-level survey depth (`depth_cont`) from the IPM photogrammetry
file; Torres Strait uses continuous logger mean depth (`depth_min_mean`) from EcoRRAP in-situ
loggers (available at Masig Reef only). Torres Strait does not have a `depth_cont` equivalent.
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
    codes = sort(unique(skipmissing(df.habitat_area)))
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

# plots — unique physical plot labels within each dataset
on_plot_n  = _nu(on_g, :plot)
ts_quad_n  = _nu(ts_g, :quadrat_number)

# depth — ON: plot-level survey depth (`depth_cont`); TS: continuous logger mean depth
# (`depth_min_mean`, available at Masig Reef only). Categorical `depth_cat` is not used here.
function _depth_stat(g, s, depth_col::Symbol)
    rows = unique(vcat(
        dropmissing(g[:, [:site_code, :habitat_area, depth_col]]),
        dropmissing(s[:, [:site_code, :habitat_area, depth_col]])
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

# habitat
on_hab_str = _hab_str(on_g)
ts_hab_str = _hab_str(ts_g)

# lon / lat (site-level)
on_lon_str = _site_stat(unique(dropmissing(on_g[:, [:site_code, :lon]])), :site_code, :lon)
on_lat_str = _site_stat(unique(dropmissing(on_g[:, [:site_code, :lat]])), :site_code, :lat)
ts_lon_str = _site_stat(unique(dropmissing(ts_g[:, [:site_code, :lon]])), :site_code, :lon)
ts_lat_str = _site_stat(unique(dropmissing(ts_g[:, [:site_code, :lat]])), :site_code, :lat)

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

# colony_id (TS photogrammetry only; ON quadrat-based so N/A)
on_cid = "N/A"
ts_cid_g = _nu(ts_g, :colony_id)
ts_cid_s = _nu(ts_s, :colony_id)

println(io, """
| Identifier | Description | Offshore North<br>Count / [Range] | Torres Strait<br>Count / [Range] |
|------------|-------------|-----------------------------------|----------------------------------|
| `diam` / `diam_mort` | Estimated coral diameter at observation (cm): `diam` for growth obs., `diam_mort` for mortality obs. | Growth: $(on_g_diam)<br>Survival: $(on_s_diam) | Growth: $(ts_g_diam)<br>Survival: $(ts_s_diam) |
| `plot` / `quadrat` | Specific monitored plot (ON) or quadrat (TS) | $(on_plot_n) | $(ts_quad_n) |
| `depth` | Continuous depth (m): ON uses plot-level survey depth (`depth_cont`) from IPM photogrammetry file; TS uses logger mean depth (`depth_min_mean`) from EcoRRAP in-situ loggers (Masig Reef only) | $(on_depth_str) | $(ts_depth_str) |
| `temp` | Mean of daily temperature (°C) from EcoRRAP in-situ loggers; daily-max mean for ON, daily mean for TS | $(on_temp_str) | $(ts_temp_str) |
| `habitat` | Reef habitat zone code (`habitat_area`) | $(on_hab_str) | $(ts_hab_str) |
| `lon` | Longitude (decimal degrees) | $(on_lon_str) | $(ts_lon_str) |
| `lat` | Latitude (decimal degrees) | $(on_lat_str) | $(ts_lat_str) |
| `site` | Reef-level site code | $(on_site_n) | $(ts_site_n) |
| `reef` | Monitored reef | $(on_reef_str) | $(ts_reef_str) |
| `taxa` | Individual taxonomic group | Growth: $(on_taxa_g_n)<br>Survival: $(on_taxa_s_n) | Growth: $(ts_taxa_g_n)<br>Survival: $(ts_taxa_s_n) |
| `functional_group` | Morphological grouping of taxa | $(on_fg_n) | $(ts_fg_n) |
| `colony_id` | Individual coral identifier | $(on_cid) | Growth: $(ts_cid_g)<br>Survival: $(ts_cid_s) |
""")

# ─── Write to file ────────────────────────────────────────────────────────────

mkpath(dirname(REPORT_FILE))
write(REPORT_FILE, String(take!(io)))
@info "Report written" file = REPORT_FILE

end # let
