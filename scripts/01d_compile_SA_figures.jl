"""
Compile exploratory sensitivity analysis figures into composite paper plots.

Requires 01a, 01b, and 01c to have been run first to produce the source PNGs.

Produces four figures:
  1. SA_growth     — (A) Offshore North, (B) Torres Strait, (C) Moore Reef, (D) Masig Reef
  2. SA_survival   — same layout for survival
  3. SA_binned_growth    — (A) Offshore North binned, (B) Torres Strait binned
  4. SA_binned_survival  — same for survival
"""

include(joinpath(@__DIR__, "common.jl"))
using FileIO

# ── Helper ─────────────────────────────────────────────────────────────────────

function load_panel(path::String)
    isfile(path) || error("Missing source figure: $path\nRun the relevant 01a/01b/01c script first.")
    return load(path)
end

function compile_composite(
    paths::Vector{String},
    labels::Vector{String};
    ncols::Int=1
)
    images = load_panel.(paths)
    nrows  = ceil(Int, length(images) / ncols)

    # Size: each column as wide as widest image in that column; each row as tall as tallest
    col_widths  = [maximum(size(images[i], 2) for i in idx:ncols:length(images)) for idx in 1:ncols]
    row_heights = [maximum(size(images[i], 1) for i in ((r-1)*ncols+1):min(r*ncols, length(images))) for r in 1:nrows]

    fig = Figure(; size=(sum(col_widths), sum(row_heights)))

    for (k, (img, lbl)) in enumerate(zip(images, labels))
        row = div(k - 1, ncols) + 1
        col = mod(k - 1, ncols) + 1
        ax  = Axis(fig[row, col]; aspect=DataAspect())
        # FileIO loads as (height × width); Makie image! expects (width × height)
        image!(ax, rotr90(img))
        hidedecorations!(ax)
        hidespines!(ax)
        Label(
            fig[row, col, TopLeft()], lbl;
            fontsize=48, font=:bold, halign=:left, valign=:top,
            padding=(6, 0, 4, 0)
        )
    end

    colgap!(fig.layout, 4)
    rowgap!(fig.layout, 4)

    return fig
end

# ── Paths ──────────────────────────────────────────────────────────────────────

sa_dir     = joinpath(FIG_DIR, "sensitivity")
output_dir = joinpath(sa_dir, "combined")
mkpath(output_dir)

on_dir    = joinpath(sa_dir, "offshore_north")
ts_dir    = joinpath(sa_dir, "torres_strait")

labels_4 = ["(A)", "(B)", "(C)", "(D)"]
labels_2 = ["(A)", "(B)"]

# ── Figure 1: Growth (regional + reef) ────────────────────────────────────────

growth_paths = [
    joinpath(on_dir, "overall", "Si_offshore_north_growth_overall.png"),       # A
    joinpath(ts_dir, "overall", "Si_torres_strait_growth_overall.png"),        # B
    joinpath(on_dir, "moore",   "Si_offshore_north_moore_growth_overall.png"), # C
    joinpath(ts_dir, "masig",   "Si_torres_strait_masig_growth_overall.png"),  # D
]

fig1 = compile_composite(growth_paths, labels_4)
display(fig1)
save(joinpath(output_dir, "SA_growth.png"), fig1; px_per_unit=DPI)
@info "Saved Figure 1: SA_growth.png"

# ── Figure 2: Survival (regional + reef) ──────────────────────────────────────

surv_paths = [
    joinpath(on_dir, "overall", "Si_offshore_north_survival_overall.png"),       # A
    joinpath(ts_dir, "overall", "Si_torres_strait_survival_overall.png"),        # B
    joinpath(on_dir, "moore",   "Si_offshore_north_moore_survival_overall.png"), # C
    joinpath(ts_dir, "masig",   "Si_torres_strait_masig_survival_overall.png"),  # D
]

fig2 = compile_composite(surv_paths, labels_4)
display(fig2)
save(joinpath(output_dir, "SA_survival.png"), fig2; px_per_unit=DPI)
@info "Saved Figure 2: SA_survival.png"

# ── Figure 3: Binned growth analysis ──────────────────────────────────────────

binned_growth_paths = [
    joinpath(on_dir, "overall", "sensitivity_growth_offshore_north.png"), # A
    joinpath(ts_dir, "overall", "sensitivity_growth_torres_strait.png"),  # B
]

fig3 = compile_composite(binned_growth_paths, labels_2)
display(fig3)
save(joinpath(output_dir, "SA_binned_growth.png"), fig3; px_per_unit=DPI)
@info "Saved Figure 3: SA_binned_growth.png"

# ── Figure 4: Binned survival analysis ────────────────────────────────────────

binned_surv_paths = [
    joinpath(on_dir, "overall", "sensitivity_survival_offshore_north.png"), # A
    joinpath(ts_dir, "overall", "sensitivity_survival_torres_strait.png"),  # B
]

fig4 = compile_composite(binned_surv_paths, labels_2)
display(fig4)
save(joinpath(output_dir, "SA_binned_survival.png"), fig4; px_per_unit=DPI)
@info "Saved Figure 4: SA_binned_survival.png"
