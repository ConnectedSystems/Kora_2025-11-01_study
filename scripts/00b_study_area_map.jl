"""
Study area map.

GBR + Torres Strait overview highlighting the EcoRRAP study reefs,
with close-up insets for each study region.

Requires the canonical reefs geopackage to have been placed in the study data
directory (already present as rrap_canonical_*.gpkg). Torres Strait study reef
locations are taken directly from `ecorrap_unified.parquet` (site_code
TSAU/TSMA/TSDU) rather than name-matched against the geopackage, since that
geopackage has no reliable entries for those three reefs.
"""

include(joinpath(@__DIR__, "common.jl"))
const GI = GDF.GeoInterface

# Requires: ]add NaturalEarth
using NaturalEarth

# ── Data paths ─────────────────────────────────────────────────────────────────
gpkg_path = joinpath(@__DIR__, "..", "data", "rrap_canonical_2025-07-15-T10-48-29.gpkg")
on_csv = joinpath(@__DIR__, "..", "data", "EcoRRAP data for IPM_250624.csv")
unified_parquet_path = joinpath(@__DIR__, "..", "data", "ecorrap_unified.parquet")

# ── Load canonical reefs ───────────────────────────────────────────────────────
reefs = GDF.read(gpkg_path)

# ── Torres Strait reef coordinates ────────────────────────────────────────────
# Name-matching the three TS study reefs against the canonical geopackage picks
# up the wrong polygons ("Akone Reef" and "Au-Masig Reef" are, per their
# coordinates, different reefs many km away from the actual survey sites; the
# geopackage has no Dungeness entry at all). Use the reef centroid recorded in
# the unified EcoRRAP dataset instead, keyed by site_code, which is authoritative.
unified_df = DataFrame(Parquet2.Dataset(unified_parquet_path))
function ts_reef_coord(site_code)
    sub = unified_df[unified_df.site_code.==site_code, :]
    isempty(sub) && error("No rows found for site_code $site_code in $unified_parquet_path")
    lons = unique(skipmissing(sub.reef_lon))
    lats = unique(skipmissing(sub.reef_lat))
    (length(lons) == 1 && length(lats) == 1) ||
        error("Ambiguous reef_lon/reef_lat for site_code $site_code")
    return (lon=lons[1], lat=lats[1])
end
const AUKANE_COORD = ts_reef_coord("TSAU")
const MASIG_COORD = ts_reef_coord("TSMA")
const DUNGENESS_COORD = ts_reef_coord("TSDU")

# ── Load EcoRRAP observation data ─────────────────────────────────────────────
on_data = CSV.read(
    on_csv, DataFrame; types=Dict(:LAT => Float64, :LONG => Float64), missingstring="NA"
)

# Unique reef names (lowercase for matching)
on_study_reefs = lowercase.(unique(on_data.Reef))

# ── Tag reefs in geopackage by role ───────────────────────────────────────────
reef_name_lc = lowercase.(reefs.reef_name)

# Normalise CSV reef names: replace underscores with spaces so "lady_musgrave"
# matches "lady musgrave reef (...)" in the geopackage.
on_study_norm = replace.(on_study_reefs, "_" => " ")

# Geographic bounds for the Offshore North study region — used to disambiguate
# reefs that share names with reefs elsewhere (e.g. "chicken" reef also exists
# in the Torres Strait).
const ON_LON_RANGE = (144.5, 148.5)
const ON_LAT_RANGE = (-18.5, -14.0)

on_in_region =
    (reefs.LON .>= ON_LON_RANGE[1]) .& (reefs.LON .<= ON_LON_RANGE[2]) .&
    (reefs.LAT .>= ON_LAT_RANGE[1]) .& (reefs.LAT .<= ON_LAT_RANGE[2])

on_reef_mask =
    [any(occursin(r, n) for r in on_study_norm) for n in reef_name_lc] .& on_in_region

# Torres Strait study reefs are NOT identified via geopackage name-matching --
# see the coordinate block above. `ts_reef_mask` stays all-false so the TS
# close-up only draws unhighlighted background reefs from the geopackage; the
# three study reefs themselves are plotted from AUKANE_COORD / MASIG_COORD /
# DUNGENESS_COORD instead.
ts_reef_mask = falses(nrow(reefs))
ts_lons = [AUKANE_COORD.lon, MASIG_COORD.lon, DUNGENESS_COORD.lon]
ts_lats = [AUKANE_COORD.lat, MASIG_COORD.lat, DUNGENESS_COORD.lat]

# The featured reef for the ON inset label (TS featured reef uses MASIG_COORD).
moore_mask = occursin.(r"(?i)moore", reefs.reef_name)

# Additional reefs requested by reviewer
lizard_mask = occursin.(r"(?i)lizard.*(island|reef)", reefs.reef_name) .& on_in_region

# ── Helper: extract exterior ring of a polygon geometry ───────────────────────
"""
    polygon_coords(geom) -> (xs, ys)

Return the exterior ring x/y vectors for a GeoInterface-compatible geometry,
or empty vectors if extraction fails (e.g. POINT, MULTIPOLYGON edge cases).
"""
function polygon_coords(geom)
    try
        gt = GI.geomtrait(geom)
        if gt isa GI.PolygonTrait
            ring = GI.getexterior(geom)
            xs = [GI.x(p) for p in GI.getpoint(ring)]
            ys = [GI.y(p) for p in GI.getpoint(ring)]
            return xs, ys
        elseif gt isa GI.MultiPolygonTrait
            # Return the largest polygon by point count
            best_xs, best_ys = Float64[], Float64[]
            for poly in GI.getgeom(geom)
                ring = GI.getexterior(poly)
                xs = [GI.x(p) for p in GI.getpoint(ring)]
                ys = [GI.y(p) for p in GI.getpoint(ring)]
                if length(xs) > length(best_xs)
                    best_xs, best_ys = xs, ys
                end
            end
            return best_xs, best_ys
        end
    catch
    end
    return Float64[], Float64[]
end

# ── Helper: draw reef polygons for a mask ─────────────────────────────────────
function draw_polygons!(ax, df, mask; color=:gray70, linewidth=0.5, fillcolor=nothing)
    for geom in df.geometry[mask]
        xs, ys = polygon_coords(geom)
        isempty(xs) && continue
        if !isnothing(fillcolor)
            poly!(ax, Point2f.(xs, ys); color=fillcolor, strokecolor=color,
                strokewidth=linewidth)
        else
            lines!(ax, xs, ys; color=color, linewidth=linewidth)
        end
    end
end

# ── Bounding boxes for close-up insets ────────────────────────────────────────
# Pad around the extent of each study group's reef centroids
function bbox_with_pad(lons, lats; pad=0.4)
    return (
        lon_min=minimum(lons) - pad,
        lon_max=maximum(lons) + pad,
        lat_min=minimum(lats) - pad,
        lat_max=maximum(lats) + pad
    )
end

@info "Matched $(sum(on_reef_mask)) ON reefs in geopackage; TS reefs plotted from unified dataset coordinates"
@info "Lizard Island matches: $(sum(lizard_mask))"

# ── Natural Earth land polygons (for close-up backgrounds + Australia inset) ──
land_50m = naturalearth("land", 50)
places_50m = naturalearth("populated_places", 50)
coast_50m = naturalearth("coastline", 50)

# Fall back to region bounds if name matching returns nothing
function safe_bbox(lons, lats, fallback_lon, fallback_lat; pad=0.4)
    isempty(lons) && return (
        lon_min=fallback_lon[1] - pad, lon_max=fallback_lon[2] + pad,
        lat_min=fallback_lat[1] - pad, lat_max=fallback_lat[2] + pad)
    return bbox_with_pad(lons, lats; pad=pad)
end

on_bbox = safe_bbox(reefs.LON[on_reef_mask], reefs.LAT[on_reef_mask],
    ON_LON_RANGE, ON_LAT_RANGE; pad=0.3)
ts_bbox = bbox_with_pad(ts_lons, ts_lats; pad=0.5)

# Expand both close-up bboxes to the same lon/lat span so panels B and C are
# plotted at the same geographic scale (same km-per-pixel with DataAspect).
let
    lon_span = max(on_bbox.lon_max - on_bbox.lon_min, ts_bbox.lon_max - ts_bbox.lon_min)
    lat_span = max(on_bbox.lat_max - on_bbox.lat_min, ts_bbox.lat_max - ts_bbox.lat_min)
    function _expand(b)
        lon_c = (b.lon_min + b.lon_max) / 2
        lat_c = (b.lat_min + b.lat_max) / 2
        return (lon_min=lon_c - lon_span / 2, lon_max=lon_c + lon_span / 2,
            lat_min=lat_c - lat_span / 2, lat_max=lat_c + lat_span / 2)
    end
    global on_bbox = _expand(on_bbox)
    global ts_bbox = _expand(ts_bbox)
end

# ── Overview extent (northern GBR + Torres Strait) ────────────────────────────
# Bounded to the latitudes/longitudes spanned by the two study regions plus a
# margin of context (Cape York to just south of Townsville); the southern GBR
# and Coral Sea carry no study reefs.
GBR_LON = (141.5, 150.0)
GBR_LAT = (-20.0, -9.0)
gbr_bbox = (lon_min=GBR_LON[1], lon_max=GBR_LON[2], lat_min=GBR_LAT[1], lat_max=GBR_LAT[2])

# ── Shared styling constants ──────────────────────────────────────────────────
COL_ON = :dodgerblue   # Offshore North
COL_TS = :darkorange   # Torres Strait
COL_ALL = (:gray75, 0.5)
STAR_SIZE = 24
STUDY_SIZE = 7
ALL_SIZE = 2

# ── Draw a north arrow in the top-right corner of an axis ────────────────────
function north_arrow!(ax; x=0.5, y=0.1, len=0.6, hw=0.012, ah=0.22, fontsize=14)
    lines!(ax, [x, x], [y, y + len - ah]; space=:relative, color=:black, linewidth=1.5)
    poly!(ax, Point2f[(x, y + len), (x - hw, y + len - ah), (x + hw, y + len - ah)];
        space=:relative, color=:black)
    return text!(ax, x, y + len + 0.03;
        text="N", space=:relative,
        align=(:center, :bottom), fontsize=fontsize, font=:bold)
end

# ── Geographic scale bar ──────────────────────────────────────────────────────
function draw_scale_bar!(ax, bbox; km=50, pad_frac=0.05, fontsize=12, halign=:left)
    lat_c = (bbox.lat_min + bbox.lat_max) / 2
    bar_deg = km / (111.0 * cosd(lat_c))
    lon_span = bbox.lon_max - bbox.lon_min
    lat_span = bbox.lat_max - bbox.lat_min
    x0 =
        halign === :right ?
        bbox.lon_max - lon_span * pad_frac - bar_deg :
        bbox.lon_min + lon_span * pad_frac
    y0 = bbox.lat_min + lat_span * pad_frac
    cap = lat_span * 0.012
    lines!(ax, [x0, x0 + bar_deg], [y0, y0]; color=:white, linewidth=5)     # halo
    lines!(ax, [x0, x0 + bar_deg], [y0, y0]; color=:black, linewidth=2)
    lines!(ax, [x0, x0], [y0 - cap, y0 + cap]; color=:black, linewidth=2)
    lines!(
        ax, [x0 + bar_deg, x0 + bar_deg], [y0 - cap, y0 + cap]; color=:black, linewidth=2
    )
    return text!(ax, x0 + bar_deg / 2, y0 + cap + 0.003;
        text="$km km", align=(:center, :bottom), fontsize=fontsize)
end

# ── Draw Natural Earth land polygons clipped to a bbox ────────────────────────
function draw_land!(ax, land_fc, bbox; strokecolor=(:black, 0.4), strokewidth=0.5)
    for feature in land_fc
        isnothing(feature.geometry) && continue
        xs, ys = polygon_coords(feature.geometry)
        isempty(xs) && continue
        any(bbox.lon_min .<= xs .<= bbox.lon_max) || continue
        any(bbox.lat_min .<= ys .<= bbox.lat_max) || continue
        poly!(ax, Point2f.(xs, ys); color=:saddlebrown,
            strokecolor=strokecolor, strokewidth=strokewidth)
    end
end

# ── City labels from Natural Earth populated_places ─────────────────────────────

# Collect coastline vertices clipped to a padded bbox.
function coast_vertices(coast_fc, bbox; pad=2.0)
    lons, lats = Float64[], Float64[]
    function _add!(geom)
        for p in GI.getpoint(geom)
            lo, la = GI.x(p), GI.y(p)
            (
                bbox.lon_min - pad <= lo <= bbox.lon_max + pad &&
                bbox.lat_min - pad <= la <= bbox.lat_max + pad
            ) || continue
            push!(lons, lo)
            push!(lats, la)
        end
    end
    for feat in coast_fc
        isnothing(feat.geometry) && continue
        gt = GI.geomtrait(feat.geometry)
        if gt isa GI.MultiLineStringTrait
            for line in GI.getgeom(feat.geometry)
                _add!(line)
            end
        else
            _add!(feat.geometry)
        end
    end
    return lons, lats
end

# Haversine distance in kilometres between two lon/lat points.
function haversine_km(lo1, la1, lo2, la2)
    R = 6371.0
    return R *
           2asin(
        sqrt(
            sin(deg2rad(la2 - la1) / 2)^2 +
            cos(deg2rad(la1)) * cos(deg2rad(la2)) * sin(deg2rad(lo2 - lo1) / 2)^2)
    )
end

"""
    draw_city_labels!(ax, places_fc, bbox; scalerank_max, fontsize, coast_fc, max_coast_km)

Plot a dot + right-aligned name label for each populated place whose SCALERANK
is ≤ `scalerank_max`, whose coordinates fall inside `bbox`, and (when `coast_fc`
is provided) whose nearest coastline vertex is within `max_coast_km` km.
"""
function draw_city_labels!(ax, places_fc, bbox; scalerank_max=7, fontsize=11,
    coast_fc=nothing, max_coast_km=100.0)
    c_lons, c_lats =
        isnothing(coast_fc) ? (Float64[], Float64[]) :
        coast_vertices(coast_fc, bbox)
    for feat in places_fc
        isnothing(feat.geometry) && continue
        lon = GI.x(feat.geometry)
        lat = GI.y(feat.geometry)
        bbox.lon_min <= lon <= bbox.lon_max || continue
        bbox.lat_min <= lat <= bbox.lat_max || continue
        feat.SCALERANK <= scalerank_max || continue
        if !isempty(c_lons)
            minimum(haversine_km(lon, lat, cl, ca)
                    for (cl, ca) in zip(c_lons, c_lats)) <= max_coast_km || continue
        end
        scatter!(ax, [lon], [lat]; color=:black, markersize=5, marker=:circle)
        text!(ax, lon - 0.15, lat;
            text=feat.NAME, fontsize=fontsize, align=(:right, :center))
    end
end

# ── Small Australia locator inset (lower-left of overview panel) ───────────────
function draw_australia_inset!(fig, position; land_fc, gbr_lon=GBR_LON, gbr_lat=GBR_LAT)
    ax_i = Axis(fig[position...];
        width=Relative(0.32), height=Relative(0.22),
        halign=0.08, valign=0.06,
        aspect=DataAspect(),
        backgroundcolor=:aliceblue)
    hidedecorations!(ax_i)
    au_bbox = (lon_min=112.0, lon_max=156.0, lat_min=-44.0, lat_max=-8.0)
    draw_land!(ax_i, land_fc, au_bbox; strokecolor=:black, strokewidth=1.2)
    bx = [gbr_lon[1], gbr_lon[2], gbr_lon[2], gbr_lon[1], gbr_lon[1]]
    by = [gbr_lat[1], gbr_lat[1], gbr_lat[2], gbr_lat[2], gbr_lat[1]]
    lines!(ax_i, bx, by; color=:red, linewidth=2)
    xlims!(ax_i, 112.0, 156.0)
    return ylims!(ax_i, -44.0, -8.0)
end

# ── Draw the overview panel ───────────────────────────────────────────────────
function draw_overview!(ax, df, on_mask, moore_mask, ts_lons, ts_lats, masig_coord)
    # All reef centroids
    scatter!(ax, df.LON, df.LAT;
        color=COL_ALL, markersize=ALL_SIZE, label="All reefs")

    # Study reef centroids
    scatter!(ax, df.LON[on_mask], df.LAT[on_mask];
        color=COL_ON, markersize=STUDY_SIZE, label="Offshore North reefs")
    scatter!(ax, ts_lons, ts_lats;
        color=COL_TS, markersize=STUDY_SIZE, label="Torres Strait reefs")

    # Featured reefs
    scatter!(ax, df.LON[moore_mask], df.LAT[moore_mask];
        color=COL_ON, marker=:star5, markersize=STAR_SIZE,
        strokecolor=:white, strokewidth=1, label="Moore Reef")
    scatter!(ax, [masig_coord.lon], [masig_coord.lat];
        color=COL_TS, marker=:star5, markersize=STAR_SIZE,
        strokecolor=:white, strokewidth=1, label="Masig Reef")

    xlims!(ax, GBR_LON...)
    return ylims!(ax, GBR_LAT...)
end

# ── Draw a close-up panel ─────────────────────────────────────────────────────
function draw_closeup!(ax, df, bbox, study_mask, featured_mask;
    study_color=COL_ON, study_label="Study reefs", feat_label="Featured reef"
)
    in_box =
        (df.LON .>= bbox.lon_min) .& (df.LON .<= bbox.lon_max) .&
        (df.LAT .>= bbox.lat_min) .& (df.LAT .<= bbox.lat_max)

    # Background reefs — try polygon outlines, scatter as fallback
    bg_mask = in_box .& .!study_mask
    if any(bg_mask)
        draw_polygons!(ax, df, bg_mask; color=(:gray60, 0.5), linewidth=0.4)
        scatter!(ax, df.LON[bg_mask], df.LAT[bg_mask];
            color=(:gray60, 0.0), markersize=0)   # invisible — just for bounds
    end

    # Study reefs
    if any(in_box .& study_mask)
        draw_polygons!(ax, df, in_box .& study_mask;
            color=study_color, fillcolor=(study_color, 0.25), linewidth=1.0)
        scatter!(ax, df.LON[in_box .& study_mask], df.LAT[in_box .& study_mask];
            color=(study_color, 0.9), markersize=9, label=study_label)
    end

    # Featured reef star
    if any(featured_mask)
        scatter!(ax, df.LON[featured_mask], df.LAT[featured_mask];
            color=study_color, marker=:star5, markersize=STAR_SIZE,
            strokecolor=:white, strokewidth=1, label=feat_label)
    end

    xlims!(ax, bbox.lon_min, bbox.lon_max)
    ylims!(ax, bbox.lat_min, bbox.lat_max)
    ax.xlabel = "Longitude"
    return ax.ylabel = "Latitude"
end

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — Overview + close-ups, no sensitive plot coordinates
# Layout: [overview (rows 1–2, col 1)] | [ON close-up (row 1, col 2)]
#                                        [TS close-up (row 2, col 2)]
# ══════════════════════════════════════════════════════════════════════════════

fig1 = Figure(; size=(1050, 900))

ax1_ov = Axis(
    fig1[1:2, 1];
    title="Study Context: Northern GBR and Torres Strait",
    xlabel="Longitude", ylabel="Latitude",
    titlesize=18, xlabelsize=18, ylabelsize=18,
    aspect=DataAspect(),
    backgroundcolor=:aliceblue
)

draw_land!(ax1_ov, land_50m, gbr_bbox)
draw_overview!(ax1_ov, reefs, on_reef_mask, moore_mask, ts_lons, ts_lats, MASIG_COORD)
if any(moore_mask)
    i = findfirst(moore_mask)
    text!(ax1_ov, reefs.LON[i], reefs.LAT[i] + 0.15;
        text="Moore", fontsize=11, align=(:center, :bottom))
end
text!(ax1_ov, MASIG_COORD.lon, MASIG_COORD.lat + 0.15;
    text="Masig", fontsize=11, align=(:center, :bottom))
draw_scale_bar!(ax1_ov, gbr_bbox; km=200, halign=:right)
draw_city_labels!(
    ax1_ov, places_50m, gbr_bbox; scalerank_max=6, coast_fc=coast_50m, max_coast_km=100.0
)
Label(fig1[1:2, 1, TopLeft()], "(A)"; fontsize=18, font=:bold, padding=(4, 0, 4, 0))
draw_australia_inset!(fig1, (1:2, 1); land_fc=land_50m)

ax1_ts = Axis(
    fig1[1, 2];
    title="Torres Strait",
    xlabel="Longitude", ylabel="Latitude",
    titlesize=18, xlabelsize=18, ylabelsize=18,
    aspect=DataAspect(),
    backgroundcolor=:aliceblue
)
draw_land!(ax1_ts, land_50m, ts_bbox)
# ts_reef_mask is all-false (see coordinate block above) so this only draws
# unhighlighted background reefs from the geopackage; the study reefs are
# plotted from the unified-dataset coordinates below.
draw_closeup!(ax1_ts, reefs, ts_bbox, ts_reef_mask, falses(nrow(reefs));
    study_color=COL_TS, study_label="Torres Strait reefs", feat_label="Masig Reef")
scatter!(ax1_ts, [AUKANE_COORD.lon, DUNGENESS_COORD.lon], [AUKANE_COORD.lat, DUNGENESS_COORD.lat];
    color=(COL_TS, 0.9), markersize=9, label="Torres Strait reefs")
text!(ax1_ts, AUKANE_COORD.lon + 0.08, AUKANE_COORD.lat;
    text="Aukane", fontsize=11, align=(:left, :center))
text!(ax1_ts, DUNGENESS_COORD.lon + 0.04, DUNGENESS_COORD.lat - 0.04;
    text="Dungeness", fontsize=11, align=(:left, :top))
scatter!(ax1_ts, [MASIG_COORD.lon], [MASIG_COORD.lat];
    color=COL_TS, marker=:star5, markersize=STAR_SIZE,
    strokecolor=:white, strokewidth=1, label="Masig Reef")
text!(ax1_ts, MASIG_COORD.lon, MASIG_COORD.lat + 0.08;
    text="Masig", fontsize=11, align=(:center, :bottom))
draw_scale_bar!(ax1_ts, ts_bbox; km=50)
Label(fig1[1, 2, TopLeft()], "(B)"; fontsize=18, font=:bold, padding=(4, 0, 4, 0))

ax1_on = Axis(
    fig1[2, 2];
    title="Offshore North",
    xlabel="Longitude", ylabel="Latitude",
    titlesize=18, xlabelsize=18, ylabelsize=18,
    aspect=DataAspect(),
    backgroundcolor=:aliceblue
)
draw_land!(ax1_on, land_50m, on_bbox)
draw_closeup!(ax1_on, reefs, on_bbox, on_reef_mask, moore_mask;
    study_color=COL_ON, study_label="Offshore North reefs", feat_label="Moore Reef")
if any(moore_mask)
    i = findfirst(moore_mask)
    text!(ax1_on, reefs.LON[i], reefs.LAT[i] + 0.06;
        text="Moore", fontsize=11, align=(:center, :bottom))
end
if any(lizard_mask)
    i = findfirst(lizard_mask)
    text!(ax1_on, reefs.LON[i], reefs.LAT[i] + 0.06;
        text="Lizard Is.", fontsize=11, align=(:center, :bottom))
end
draw_scale_bar!(ax1_on, on_bbox; km=50)
Label(fig1[2, 2, TopLeft()], "(C)"; fontsize=18, font=:bold, padding=(4, 0, 4, 0))

colsize!(fig1.layout, 1, Auto(0.5))
colsize!(fig1.layout, 2, Auto(0.25))

Legend(
    fig1[3, 1],
    [
        MarkerElement(; color=COL_ON, marker=:circle, markersize=STUDY_SIZE),
        MarkerElement(; color=COL_TS, marker=:circle, markersize=STUDY_SIZE),
        MarkerElement(; color=COL_ON, marker=:star5, markersize=STAR_SIZE,
            strokecolor=:white, strokewidth=1),
        MarkerElement(; color=COL_TS, marker=:star5, markersize=STAR_SIZE,
            strokecolor=:white, strokewidth=1)
    ],
    ["Offshore North reefs", "Torres Strait reefs", "Moore Reef", "Masig Reef"];
    orientation=:horizontal, framevisible=false, nbanks=2, labelsize=14
)

ax1_north = Axis(fig1[3, 2])
hidedecorations!(ax1_north)
hidespines!(ax1_north)
xlims!(ax1_north, 0, 1)
ylims!(ax1_north, 0, 1)
north_arrow!(ax1_north; fontsize=16)

fig1_path = joinpath(FIG_DIR, "study_area_overview.png")
save(fig1_path, fig1; px_per_unit=DPI)
@info "Saved Figure 1: $fig1_path"
