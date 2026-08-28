
include(joinpath(@__DIR__, "common.jl"))

# ── Paths to cached PAWN results ──────────────────────────────────────────────
moore_data_dir = joinpath(
    OUTPUT_DIR, "sensitivity", "offshore_north", "moore", "ensemble"
)
masig_data_dir = joinpath(
    OUTPUT_DIR, "sensitivity", "torres_strait", "masig", "ensemble"
)

moore_cons_pawn  = load_result(joinpath(moore_data_dir,  "16071S_constrained_pawn_results.h5"))
moore_unc_pawn   = load_result(joinpath(moore_data_dir,  "16071S_unconstrained_pawn_results.h5"))
masig_cons_pawn  = load_result(joinpath(masig_data_dir,  "masig_constrained_pawn_results.h5"))
masig_unc_pawn   = load_result(joinpath(masig_data_dir,  "masig_unconstrained_pawn_results.h5"))

combined_fig_dir = joinpath(FIG_DIR, "sensitivity", "combined")
mkpath(combined_fig_dir)

# ── Helper: plot PAWN heatmap into a GridLayout sub-position ──────────────────
"""
    plot_pawn_heatmap_into!(layout_pos, Si, title; ...)

Draw a PAWN sensitivity heatmap into `layout_pos` (a Makie GridLayout position
or sub-layout cell).  Returns the `HeatmapPlot` object so the caller can attach
a shared colorbar.

`colorrange` should be passed explicitly when panels are combined in one figure
so that both heatmaps share a consistent scale.
"""
function plot_pawn_heatmap_into!(
    layout_pos,
    Si::AbstractDimArray,
    title::String;
    stats::Vector{Symbol}=[:mean, :std],
    xticklabelrotation::Real=π / 2,
    colorrange::Tuple{<:Real,<:Real}=(-0.1, 0.5)
)
    factor_order  = sortperm(collect(Si[PAWNᵢ=At(:mean)]); rev=true)
    factor_labels = string.(collect(dims(Si, 1)))[factor_order]

    data       = hcat([collect(Si[PAWNᵢ=At(s)])[factor_order] for s in stats]...)
    stat_labels = string.(stats)

    ax = Axis(layout_pos[1, 1])
    hm = heatmap!(ax, data; colorrange=colorrange, colormap=:viridis)

    ax.xticks             = (1:length(factor_labels), factor_labels)
    ax.yticks             = (1:length(stat_labels),   stat_labels)
    ax.yreversed          = true
    ax.xticklabelrotation = xticklabelrotation
    ax.xlabelsize         = 14
    ax.ylabelsize         = 14
    ax.xticklabelsize     = 11
    ax.yticklabelsize     = 12
    ax.titlesize          = 14
    ax.title              = title

    return hm
end

# ── Helper: compute max PAWN value across a set of results ────────────────────
function pawn_data_max(Si::AbstractDimArray, stats::Vector{Symbol}=[:mean, :std])
    data = hcat([collect(Si[PAWNᵢ=At(s)]) for s in stats]...)
    return maximum(data)
end

# ── Publication theme ─────────────────────────────────────────────────────────
set_theme!(Theme(; fontsize=14))

stats = [:mean, :std]

# ═══════════════════════════════════════════════════════════════════════════════
# Figure A — Constrained SA (2×1: Moore top, Masig bottom)
# ═══════════════════════════════════════════════════════════════════════════════

cons_max = max(
    pawn_data_max(moore_cons_pawn, stats),
    pawn_data_max(masig_cons_pawn, stats),
    0.1
)
cons_cr = (-0.1, cons_max)

fig_cons = Figure(; size=(900, 560))

gl_moore_cons = fig_cons[1, 1] = GridLayout()
gl_masig_cons = fig_cons[2, 1] = GridLayout()

hm_cons = plot_pawn_heatmap_into!(
    gl_moore_cons, moore_cons_pawn, "Constrained SA - Moore Reef";
    stats=stats, colorrange=cons_cr
)
plot_pawn_heatmap_into!(
    gl_masig_cons, masig_cons_pawn, "Constrained SA - Masig";
    stats=stats, colorrange=cons_cr
)

Colorbar(fig_cons[1:2, 2], hm_cons; label="PAWN Index", width=14)

save(joinpath(combined_fig_dir, "combined_constrained_sa.png"), fig_cons; px_per_unit=DPI)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure B — Unconstrained SA (2×1: Moore top, Masig bottom)
# ═══════════════════════════════════════════════════════════════════════════════

unc_max = max(
    pawn_data_max(moore_unc_pawn, stats),
    pawn_data_max(masig_unc_pawn, stats),
    0.1
)
unc_cr = (-0.1, unc_max)

fig_unc = Figure(; size=(900, 560))

gl_moore_unc = fig_unc[1, 1] = GridLayout()
gl_masig_unc = fig_unc[2, 1] = GridLayout()

hm_unc = plot_pawn_heatmap_into!(
    gl_moore_unc, moore_unc_pawn, "Unconstrained SA - Moore Reef";
    stats=stats, colorrange=unc_cr
)
plot_pawn_heatmap_into!(
    gl_masig_unc, masig_unc_pawn, "Unconstrained SA - Masig";
    stats=stats, colorrange=unc_cr
)

Colorbar(fig_unc[1:2, 2], hm_unc; label="PAWN Index", width=14)

save(joinpath(combined_fig_dir, "combined_unconstrained_sa.png"), fig_unc; px_per_unit=DPI)
