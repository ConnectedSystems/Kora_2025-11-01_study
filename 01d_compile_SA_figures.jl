"""
Compile exploratory sensitivity analysis figures into composite paper plots.

Requires 01a, 01b, and 01c to have been run first to produce the source PNGs.

Produces five figures:
  1. diameter_bins_comparison     — (A) Offshore North binned, (B) Torres Strait binned
  2. SA_offshore_north_growth     — (A) Offshore North overall, (B) Moore Reef
  3. SA_offshore_north_survival   — same for survival
  4. SA_torres_strait_growth      — (A) Torres Strait overall, (B) Masig Reef
  5. SA_torres_strait_survival    — same for survival
"""

include("common.jl")
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

sa_dir = joinpath(FIG_DIR, "sensitivity")

region_configs = [
    (region="offshore_north", reef_id="moore",  reef_name="Moore"),
    (region="torres_strait",  reef_id="masig",  reef_name="Masig"),
]

output_dir = joinpath(sa_dir, "combined")
mkpath(output_dir)

# ── Figure 1: Diameter-bins comparison (Offshore North vs Torres Strait) ───────

bin_paths = [
    joinpath(sa_dir, "offshore_north", "overall", "sensitivity_growth_offshore_north.png"),
    joinpath(sa_dir, "torres_strait",  "overall", "sensitivity_growth_torres_strait.png"),
]
bin_labels = ["(A)", "(B)"]

fig_bins_growth = compile_composite(bin_paths, bin_labels)
display(fig_bins_growth)
save(joinpath(output_dir, "diameter_bins_growth_comparison.png"), fig_bins_growth; px_per_unit=DPI)

bin_surv_paths = [
    joinpath(sa_dir, "offshore_north", "overall", "sensitivity_survival_offshore_north.png"),
    joinpath(sa_dir, "torres_strait",  "overall", "sensitivity_survival_torres_strait.png"),
]

fig_bins_surv = compile_composite(bin_surv_paths, bin_labels)
display(fig_bins_surv)
save(joinpath(output_dir, "diameter_bins_survival_comparison.png"), fig_bins_surv; px_per_unit=DPI)

@info "Saved diameter-bins comparison figures"

# ── Figures 2–5: Region + reef composites (growth and survival) ────────────────

for (; region, reef_id, reef_name) in region_configs
    region_dir = joinpath(sa_dir, region)
    scale_fn   = "$(region)_$(reef_id)"

    growth_paths = [
        joinpath(region_dir, "overall", "Si_$(region)_growth_overall.png"),    # (A) region overall
        joinpath(region_dir, reef_id,   "Si_$(scale_fn)_growth_overall.png"),  # (B) reef-specific
    ]

    surv_paths = [
        joinpath(region_dir, "overall", "Si_$(region)_survival_overall.png"),
        joinpath(region_dir, reef_id,   "Si_$(scale_fn)_survival_overall.png"),
    ]

    labels = ["(A)", "(B)"]

    fig_growth = compile_composite(growth_paths, labels)
    fig_surv   = compile_composite(surv_paths,   labels)

    display(fig_growth)
    display(fig_surv)

    save(joinpath(output_dir, "SA_$(scale_fn)_growth_composite.png"),   fig_growth; px_per_unit=DPI)
    save(joinpath(output_dir, "SA_$(scale_fn)_survival_composite.png"), fig_surv;   px_per_unit=DPI)

    @info "Saved growth/survival composites for $(reef_name) ($(region))"
end
