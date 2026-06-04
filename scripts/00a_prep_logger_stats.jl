"""
Process all EcoRRAP oceanographic raw data files into period-level summary statistics
per site and depth category.

Output: data/ecorrap_logger/ocn_annual_stats.parquet

Each row represents one (site_code, depth_category, period_year) combination.
Statistics are computed as both mean and median over the chosen temporal window.

window_type options:
  :survey_year   — May–April (default); label = May start year; aligns with inter-survey period
  :water_year    — Nov–Oct; label = Oct end year
  :calendar_year — Jan–Dec; label = calendar year
  :wet_season    — Nov–Apr only; label = year of the January
  :dry_season    — May–Oct only; label = calendar year
"""

using CSV
using DataFrames
using Parquet2
using Statistics
using Dates

const OCN_DIR = joinpath(@__DIR__, "..", "data", "ecorrap_logger")
const OUTPUT_FILE = joinpath(OCN_DIR, "ocn_annual_stats.parquet")
const WINDOW_TYPE = :survey_year

# Minimum days of data required in a period before emitting a warning
const MIN_COVERAGE_DAYS = 90

# ─── Period assignment ────────────────────────────────────────────────────────

"""Return the survey year (May–April) label for a date. May 2021–Apr 2022 → 2021."""
survey_year(d::Date)::Int = month(d) >= 5 ? year(d) : year(d) - 1

"""Return the water year (Nov–Oct) label for a date. Nov 2021–Oct 2022 → 2022."""
water_year(d::Date)::Int = month(d) >= 11 ? year(d) + 1 : year(d)

"""
Assign a period label to a date under the chosen window_type.
Returns `missing` if the date falls outside the active window
(only relevant for :wet_season and :dry_season).
"""
function assign_period(d::Date, window::Symbol)::Union{Int,Missing}
    if window === :survey_year
        return survey_year(d)
    elseif window === :water_year
        return water_year(d)
    elseif window === :calendar_year
        return year(d)
    elseif window === :wet_season
        month(d) ∈ (11, 12, 1, 2, 3, 4) || return missing
        return month(d) >= 11 ? year(d) + 1 : year(d)
    elseif window === :dry_season
        month(d) ∈ (5, 6, 7, 8, 9, 10) || return missing
        return year(d)
    else
        error(
            "Unknown window_type: $window. Choose from :survey_year, :water_year, " *
            ":calendar_year, :wet_season, :dry_season"
        )
    end
end

"""Classify instrument nominal depth as shallow (S) or deep (D)."""
depth_category(d::Real)::String = d >= 8.0 ? "D" : "S"

# ─── File I/O helpers ─────────────────────────────────────────────────────────

"""
Parse the `# key : value` metadata header from an oceanographic CSV file.
Stops at the first non-comment line.
"""
function parse_ocn_header(filepath::String)::Dict{String,String}
    meta = Dict{String,String}()
    open(filepath) do fh
        for line in eachline(fh)
            startswith(line, "#") || break
            m = match(r"^#\s*([^:]+?)\s*:\s*(.+)$", line)
            isnothing(m) || (meta[strip(m[1])] = strip(m[2]))
        end
    end
    return meta
end

"""
Read the CSV data portion of an oceanographic file, skipping all `#` header lines.
Returns an empty DataFrame if no data rows are found.
"""
function read_ocn_csv(filepath::String)::DataFrame
    lines = readlines(filepath)
    data_start = findfirst(!startswith("#"), lines)
    isnothing(data_start) && return DataFrame()
    return CSV.read(
        IOBuffer(join(lines[data_start:end], "\n")),
        DataFrame;
        missingstring=["", "nan", "NaN", "NaT", "-9999", "-9999.0"]
    )
end

"""
Extract (site_code, depth_category, lat, lon) from a parsed header Dict.
site_code is the first element of platform_code before the first underscore.
"""
function extract_header_meta(meta::Dict{String,String})
    platform = get(meta, "platform_code", "")
    site_code = isempty(platform) ? "" : split(platform, "_")[1]
    depth_m = tryparse(Float64, get(meta, "instrument_nominal_depth", ""))
    lat = tryparse(Float64, get(meta, "geospatial_lat_max", ""))
    lon = tryparse(Float64, get(meta, "geospatial_lon_max", ""))
    depth_cat = isnothing(depth_m) ? "" : depth_category(depth_m)
    return (site_code=site_code, depth_cat=depth_cat, lat=lat, lon=lon)
end

"""
Recursively find all CSV files under `base_dir`, optionally filtered by `file_filter(filename)`.
"""
function find_csv_files(base_dir::String, file_filter=nothing)::Vector{String}
    files = String[]
    isdir(base_dir) || return files
    for (root, _, filenames) in walkdir(base_dir)
        for fn in filenames
            endswith(fn, ".csv") || continue
            (isnothing(file_filter) || file_filter(fn)) && push!(files, joinpath(root, fn))
        end
    end
    return files
end

# ─── Two-step aggregation helpers ─────────────────────────────────────────────
#
# Step 1: average within (site_code, depth_cat, date) across multiple instruments
#         so no site is over-weighted by number of deployed loggers.
# Step 2: aggregate those daily values to (site_code, depth_cat, period_year) statistics.

"""
Given a long DataFrame with columns [:site_code, :depth_cat, :date, value_cols...],
return daily averages within each (site, depth, date) group.
"""
function to_daily(df::DataFrame, value_cols::Vector{Symbol})::DataFrame
    isempty(df) && return df
    gdf = groupby(df, [:site_code, :depth_cat, :date])
    agg_pairs = [col => mean ∘ skipmissing => col for col in value_cols]
    return combine(gdf, agg_pairs...)
end

"""
Aggregate a daily DataFrame to period statistics.
Returns a DataFrame with columns:
  site_code, depth_cat, period_year, <col>_mean, <col>_median for each value column,
  plus n_days (coverage count).
"""
function to_period_stats(
    daily_df::DataFrame,
    value_cols::Vector{Symbol},
    window::Symbol
)::DataFrame
    isempty(daily_df) && return DataFrame()

    daily_df = copy(daily_df)
    daily_df.period_year = [assign_period(d, window) for d in daily_df.date]
    filter!(:period_year => !ismissing, daily_df)
    daily_df.period_year = Int.(daily_df.period_year)

    gdf = groupby(daily_df, [:site_code, :depth_cat, :period_year])

    agg_pairs = Pair[]
    for col in value_cols
        push!(agg_pairs, col => (x -> (v = collect(skipmissing(x)); isempty(v) ? missing : mean(v))) => Symbol("$(col)_mean"))
        push!(agg_pairs, col => (x -> (v = collect(skipmissing(x)); isempty(v) ? missing : median(v))) => Symbol("$(col)_median"))
    end
    push!(agg_pairs, first(value_cols) => length => :n_days)

    result = combine(gdf, agg_pairs...)

    # Warn on sparse coverage
    for row in eachrow(result)
        if row.n_days < MIN_COVERAGE_DAYS
            @warn "Sparse coverage" site = row.site_code depth = row.depth_cat period =
                row.period_year n_days = row.n_days min_expected = MIN_COVERAGE_DAYS
        end
    end

    return result
end

# ─── Variable-specific loaders ────────────────────────────────────────────────

"""
Load and aggregate daily TEMP_STATS files.
Returns columns: site_code, depth_cat, period_year,
  temp_mean_mean, temp_mean_median, temp_max_mean, temp_max_median, n_days_temp
"""
function process_temp_stats(window::Symbol)::DataFrame
    subdir = joinpath(OCN_DIR, "TEMP_STATS")
    files = find_csv_files(subdir)
    @info "TEMP_STATS: found $(length(files)) files"

    frames = DataFrame[]
    for fp in files
        try
            meta = parse_ocn_header(fp)
            info = extract_header_meta(meta)
            isempty(info.site_code) && continue

            df = read_ocn_csv(fp)
            isempty(df) && continue
            :MEAN ∈ propertynames(df) && :MAX ∈ propertynames(df) || continue

            df.TIME = DateTime.(df.TIME)
            df.date = Date.(df.TIME)
            df.site_code = fill(info.site_code, nrow(df))
            df.depth_cat = fill(info.depth_cat, nrow(df))

            push!(frames, select(df, [:site_code, :depth_cat, :date, :MEAN, :MAX]))
        catch e
            @warn "Skipping TEMP file" file = basename(fp) error = e
        end
    end

    isempty(frames) && return DataFrame()
    combined = reduce(vcat, frames)
    rename!(combined, :MEAN => :temp_mean, :MAX => :temp_max)

    daily = to_daily(combined, [:temp_mean, :temp_max])
    stats = to_period_stats(daily, [:temp_mean, :temp_max], window)
    rename!(stats, :n_days => :n_days_temp)
    return stats
end

"""
Load and aggregate daily PSAL_STATS files.
Returns columns: site_code, depth_cat, period_year,
  psal_mean_mean, psal_mean_median, n_days_psal
"""
function process_psal_stats(window::Symbol)::DataFrame
    subdir = joinpath(OCN_DIR, "PSAL_STATS")
    files = find_csv_files(subdir)
    @info "PSAL_STATS: found $(length(files)) files"

    frames = DataFrame[]
    for fp in files
        try
            meta = parse_ocn_header(fp)
            info = extract_header_meta(meta)
            isempty(info.site_code) && continue

            df = read_ocn_csv(fp)
            isempty(df) && continue
            :MEAN ∈ propertynames(df) || continue

            df.TIME = DateTime.(df.TIME)
            df.date = Date.(df.TIME)
            df.site_code = fill(info.site_code, nrow(df))
            df.depth_cat = fill(info.depth_cat, nrow(df))

            push!(frames, select(df, [:site_code, :depth_cat, :date, :MEAN]))
        catch e
            @warn "Skipping PSAL file" file = basename(fp) error = e
        end
    end

    isempty(frames) && return DataFrame()
    combined = reduce(vcat, frames)
    rename!(combined, :MEAN => :psal_mean)

    daily = to_daily(combined, [:psal_mean])
    stats = to_period_stats(daily, [:psal_mean], window)
    rename!(stats, :n_days => :n_days_psal)
    return stats
end

"""
Load and aggregate hourly CSPD_STATS files.
Hourly observations are averaged to daily before period aggregation.
Returns columns: site_code, depth_cat, period_year,
  cspd_mean_mean, cspd_mean_median, n_days_cspd
"""
function process_cspd_stats(window::Symbol)::DataFrame
    subdir = joinpath(OCN_DIR, "CSPD_STATS")
    files = find_csv_files(subdir)
    @info "CSPD_STATS: found $(length(files)) files"

    frames = DataFrame[]
    for fp in files
        try
            meta = parse_ocn_header(fp)
            info = extract_header_meta(meta)
            isempty(info.site_code) && continue

            df = read_ocn_csv(fp)
            isempty(df) && continue
            :MEAN ∈ propertynames(df) || continue

            df.TIME = DateTime.(df.TIME)
            df.date = Date.(df.TIME)
            df.site_code = fill(info.site_code, nrow(df))
            df.depth_cat = fill(info.depth_cat, nrow(df))

            push!(frames, select(df, [:site_code, :depth_cat, :date, :MEAN]))
        catch e
            @warn "Skipping CSPD file" file = basename(fp) error = e
        end
    end

    isempty(frames) && return DataFrame()
    combined = reduce(vcat, frames)
    rename!(combined, :MEAN => :cspd_mean)

    daily = to_daily(combined, [:cspd_mean])
    stats = to_period_stats(daily, [:cspd_mean], window)
    rename!(stats, :n_days => :n_days_cspd)
    return stats
end

"""
Load and aggregate WAVES_CSV files (W-prefix only; X-prefix files are excluded).
Sub-hourly observations are averaged to daily before period aggregation.
Uses WSSH (spectral significant wave height).
Returns columns: site_code, depth_cat, period_year,
  wave_hs_mean_mean, wave_hs_mean_median, n_days_waves
"""
function process_waves(window::Symbol)::DataFrame
    subdir = joinpath(OCN_DIR, "WAVES_CSV")
    # Include only standard W-prefix wave files; skip X-prefix variants
    files = find_csv_files(subdir, fn -> contains(fn, "_W_") && endswith(fn, "_waves.csv"))
    @info "WAVES_CSV (W only): found $(length(files)) files"

    frames = DataFrame[]
    for fp in files
        try
            meta = parse_ocn_header(fp)
            info = extract_header_meta(meta)
            isempty(info.site_code) && continue

            df = read_ocn_csv(fp)
            isempty(df) && continue
            :WSSH ∈ propertynames(df) || continue

            df.TIME = DateTime.(df.TIME)
            df.date = Date.(df.TIME)
            df.site_code = fill(info.site_code, nrow(df))
            df.depth_cat = fill(info.depth_cat, nrow(df))

            push!(frames, select(df, [:site_code, :depth_cat, :date, :WSSH]))
        catch e
            @warn "Skipping WAVES file" file = basename(fp) error = e
        end
    end

    isempty(frames) && return DataFrame()
    combined = reduce(vcat, frames)
    rename!(combined, :WSSH => :wave_hs)

    daily = to_daily(combined, [:wave_hs])
    stats = to_period_stats(daily, [:wave_hs], window)
    rename!(stats, :n_days => :n_days_waves)
    return stats
end

"""
Load and aggregate PAR daily DLI files (*_dli.csv suffix only).
Uses dli_sum_deltat (full-day Daily Light Integral, mol m⁻² day⁻¹).
Returns columns: site_code, depth_cat, period_year,
  par_dli_mean, par_dli_median, n_days_par
"""
function process_par_dli(window::Symbol)::DataFrame
    subdir = joinpath(OCN_DIR, "PAR_CSV")
    files = find_csv_files(subdir, fn -> endswith(fn, "_dli.csv"))
    @info "PAR_CSV (dli only): found $(length(files)) files"

    frames = DataFrame[]
    for fp in files
        try
            meta = parse_ocn_header(fp)
            info = extract_header_meta(meta)
            isempty(info.site_code) && continue

            df = read_ocn_csv(fp)
            isempty(df) && continue
            :dli_sum_deltat ∈ propertynames(df) || continue

            df.TIME = DateTime.(df.TIME)
            df.date = Date.(df.TIME)
            df.site_code = fill(info.site_code, nrow(df))
            df.depth_cat = fill(info.depth_cat, nrow(df))

            push!(frames, select(df, [:site_code, :depth_cat, :date, :dli_sum_deltat]))
        catch e
            @warn "Skipping PAR file" file = basename(fp) error = e
        end
    end

    isempty(frames) && return DataFrame()
    combined = reduce(vcat, frames)
    rename!(combined, :dli_sum_deltat => :par_dli)

    # DLI of 0 on the first/last partial deployment day is a recording artefact; keep
    # genuine zeros (overcast days are valid) but drop missing.
    daily = to_daily(combined, [:par_dli])
    stats = to_period_stats(daily, [:par_dli], window)
    rename!(stats, :n_days => :n_days_par)
    return stats
end

"""
Load and aggregate daily DEPTH_STATS files.
Returns columns: site_code, depth_cat, period_year,
  depth_min_mean, depth_min_median, depth_max_mean, depth_max_median,
  depth_range_mean, depth_range_median, n_days_depth
"""
function process_depth_stats(window::Symbol)::DataFrame
    subdir = joinpath(OCN_DIR, "DEPTH_STATS")
    files = find_csv_files(subdir)
    @info "DEPTH_STATS: found $(length(files)) files"

    frames = DataFrame[]
    for fp in files
        try
            meta = parse_ocn_header(fp)
            info = extract_header_meta(meta)
            isempty(info.site_code) && continue

            df = read_ocn_csv(fp)
            isempty(df) && continue
            :MIN ∈ propertynames(df) && :MAX ∈ propertynames(df) && :RANGE ∈ propertynames(df) || continue

            df.TIME = DateTime.(df.TIME)
            df.date = Date.(df.TIME)
            df.site_code = fill(info.site_code, nrow(df))
            df.depth_cat = fill(info.depth_cat, nrow(df))

            push!(frames, select(df, [:site_code, :depth_cat, :date, :MIN, :MAX, :RANGE]))
        catch e
            @warn "Skipping DEPTH file" file = basename(fp) error = e
        end
    end

    isempty(frames) && return DataFrame()
    combined = reduce(vcat, frames)
    rename!(combined, :MIN => :depth_min, :MAX => :depth_max, :RANGE => :depth_range)

    daily = to_daily(combined, [:depth_min, :depth_max, :depth_range])
    stats = to_period_stats(daily, [:depth_min, :depth_max, :depth_range], window)
    rename!(stats, :n_days => :n_days_depth)
    return stats
end

# ─── Main ─────────────────────────────────────────────────────────────────────

@info "Processing oceanographic data with window_type = :$(WINDOW_TYPE)"

temp_stats = process_temp_stats(WINDOW_TYPE)
psal_stats = process_psal_stats(WINDOW_TYPE)
cspd_stats = process_cspd_stats(WINDOW_TYPE)
wave_stats = process_waves(WINDOW_TYPE)
par_stats = process_par_dli(WINDOW_TYPE)
depth_stats = process_depth_stats(WINDOW_TYPE)

# Full outer join on (site_code, depth_cat, period_year) so every site appears
# even when some variables are absent
join_key = [:site_code, :depth_cat, :period_year]

all_stats = foldl(
    (a, b) -> outerjoin(a, b; on=join_key),
    [temp_stats, psal_stats, cspd_stats, wave_stats, par_stats, depth_stats]
)

# Record which window type produced this file
all_stats.window_type = fill(string(WINDOW_TYPE), nrow(all_stats))

# Sort for readability
sort!(all_stats, [:site_code, :depth_cat, :period_year])

@info "Writing oceanographic stats" rows = nrow(all_stats) file = OUTPUT_FILE
Parquet2.writefile(OUTPUT_FILE, all_stats)

@info "Done. Output: $OUTPUT_FILE"
