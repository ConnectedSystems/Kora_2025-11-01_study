"""
Sensitivity analysis across factor space and regions
"""

include(joinpath(@__DIR__, "common.jl"))

Random.seed!(42)

"""
    binned_sa(
        data::DataFrame,
        y_values::Vector,
        diameter_col::Symbol,
        n_bins::Int=10
    )

Perform PAWN sensitivity analysis across diameter bins.

Returns (g_bin, bin_details, per_bin_sample) where:
- g_bin: Matrix of sensitivity indices (n_bins × n_features)
- bin_details: Matrix of bin information (n_bins × 3) containing [start, mean, end]
- per_bin_sample: Minimum number of samples per bin
"""
function binned_sa(
    data::DataFrame,
    y_values::Vector,
    diameters::Vector;
    n_bins::Int=10
)
    n_obs = nrow(data)
    per_bin_sample = n_obs ÷ n_bins
    bin_ids = Kora.adaptive_min_sample_binning(diameters, per_bin_sample)

    # Drop columns that are constant across the full dataset — pawn's quantile
    # step will fail with an empty data vector if all values in a bin are identical.
    const_cols = [n for n in names(data) if length(unique(data[!, n])) == 1]
    if !isempty(const_cols)
        @warn "Dropping constant columns before sensitivity analysis" const_cols
        select!(data, Not(const_cols))
    end

    n_bins_actual = length(unique(bin_ids))
    g_bin = Matrix{Float64}(undef, n_bins_actual, ncol(data))
    bin_details = Matrix{Float64}(undef, n_bins_actual, 3)

    for i in sort(unique(bin_ids))
        bin_sel = bin_ids .== i

        bin_mean = mean(diameters[bin_sel])
        (bin_start, bin_end) = extrema(diameters[bin_sel])
        bin_details[i, :] .= bin_start, bin_mean, bin_end

        Si = pawn(data[bin_sel, :], y_values[bin_sel]; S=10)[PAWNᵢ=At(:mean)]

        g_bin[i, :] = Si
    end

    return g_bin, bin_details, per_bin_sample
end

"""
    plot_sensitivity_heatmap(
        g_bin::Matrix,
        bin_details::Matrix,
        feature_names::Vector{String},
        region::String,
        analysis_type::String,
        per_bin_sample::Int
    )

Create sensitivity analysis heatmap.
"""
function plot_sensitivity_heatmap(
    g_bin::Matrix,
    bin_details::Matrix,
    feature_names::Vector{String},
    region::String,
    analysis_type::String,
    per_bin_sample::Int
)
    n_bins = size(g_bin, 1)

    # Sort features by mean influence across bins, most influential at top
    order = sortperm(vec(mean(g_bin; dims=1)))
    g_bin_sorted = g_bin[:, order]
    feature_names_sorted = feature_names[order]

    f = Figure(; size=(800, 600))
    ax = Axis(f[1, 1])

    heatmap!(ax, g_bin_sorted)

    ax.yticks = (1:length(feature_names_sorted), feature_names_sorted)
    ax.xticks = (1:n_bins, string.(round.(bin_details[:, 2]; digits=2)))
    ax.title = "$(titlecase(replace(region, "_" => " "))) - $(analysis_type)"
    ax.ylabel = "Factors"
    ax.xlabel = "Mean Diameter of Bin\n($(per_bin_sample) samples per bin)"

    Colorbar(f[1, 2]; limits=(-0.1, max(maximum(g_bin_sorted), 1.0)), label="PAWN Index")

    return f
end

"""
    prepare_growth_data(model_results)

Prepare growth data for sensitivity analysis.
"""
function prepare_growth_data(model_results)
    all_growth = vcat(values(model_results.growth_groupings)...)
    all_y_growth = all_growth.est_1yo_growth
    diameters = Float64.(all_growth.diam)

    ignore_cols = [g for g in growth_ignore_cols if g in propertynames(all_growth)]
    select!(all_growth, Not(ignore_cols))
    cleanup_features!(all_growth)
    rename_for_display!(all_growth)

    return all_growth, all_y_growth, diameters
end

"""
    prepare_survival_data(model_results)

Prepare survival data for sensitivity analysis.
"""
function prepare_survival_data(model_results)
    all_surv = vcat(values(model_results.survival_groupings)...)
    all_y_surv = all_surv.surv
    all_y_surv[ismissing.(all_y_surv)] .= 0
    all_y_surv = Int64.(all_y_surv)
    diameters = Float64.(all_surv.diam_mort)

    ignore_cols = [g for g in surv_ignore_cols if g in propertynames(all_surv)]
    select!(all_surv, Not(ignore_cols))
    cleanup_features!(all_surv)
    rename_for_display!(all_surv)

    return all_surv, all_y_surv, diameters
end

"""
    process_region_sensitivity(
        ecorrap_file::String,
        species_file::String,
        region::String;
        n_bins::Int=10
    )

Process sensitivity analysis for growth and survival in a single region.
"""
function process_region_sensitivity(
    ecorrap_file::String,
    species_file::String,
    region::String;
    n_bins::Int=10
)
    @info "Processing region: $region"

    # Fit models
    model_results = Kora.process_ecorrap_models(
        ecorrap_file,
        species_file;
        region=region,
        save_models=false,
        plot_validation=false,
        growth_degree=2
    )

    # Prepare data
    @info "Preparing growth data..."
    all_growth, all_y_growth, growth_diams = prepare_growth_data(model_results)

    @info "Preparing survival data..."
    all_surv, all_y_surv, surv_diams = prepare_survival_data(model_results)

    # Replace any remaining NaNs across all float columns with a sentinel value.
    # NaNs arise from cleanup_features! (missing → NaN) for sparsely covered columns;
    # pawn's quantile step cannot handle them.
    for col in names(all_growth)
        eltype(all_growth[!, col]) <: AbstractFloat && replace!(all_growth[!, col], NaN => -1.0)
    end
    for col in names(all_surv)
        eltype(all_surv[!, col]) <: AbstractFloat && replace!(all_surv[!, col], NaN => -1.0)
    end

    # Perform sensitivity analyses
    @info "Analyzing growth sensitivity..."
    g_bin_growth, bin_details_growth, per_bin_sample_growth = binned_sa(
        all_growth, all_y_growth, growth_diams; n_bins=n_bins
    )

    @info "Analyzing survival sensitivity..."
    g_bin_surv, bin_details_surv, per_bin_sample_surv = binned_sa(
        all_surv, all_y_surv, surv_diams; n_bins=n_bins
    )

    # Create plots
    fig_growth = plot_sensitivity_heatmap(
        g_bin_growth, bin_details_growth, names(all_growth),
        region, "Growth Sensitivity", per_bin_sample_growth
    )

    fig_surv = plot_sensitivity_heatmap(
        g_bin_surv, bin_details_surv, names(all_surv),
        region, "Survival Sensitivity", per_bin_sample_surv
    )

    return (
        growth=(
            indices=g_bin_growth,
            bins=bin_details_growth,
            features=names(all_growth),
            figure=fig_growth
        ),
        survival=(
            indices=g_bin_surv,
            bins=bin_details_surv,
            features=names(all_surv),
            figure=fig_surv
        )
    )
end

# Main analysis
# Each region uses its own EcoRRAP data file
region_data = [
    ("offshore_north", joinpath(OUTPUT_DIR, "ecorrap_unified.parquet")),
    ("torres_strait", joinpath(OUTPUT_DIR, "ecorrap_unified.parquet"))
]
species_file = joinpath(OUTPUT_DIR, "ecorrap_to_cscape_species.csv")

results = Dict{String,NamedTuple}()

for (region, ecorrap_file) in region_data
    results[region] = process_region_sensitivity(
        ecorrap_file, species_file, region; n_bins=10
    )
end

# Display figures
for (region, _) in region_data
    @info "Displaying results for $region"
    display(results[region].growth.figure)
    display(results[region].survival.figure)
end

# Save results
for (region, _) in region_data
    region_overall_dir = joinpath(FIG_DIR, "sensitivity", region, "overall")
    mkpath(region_overall_dir)
    save(
        joinpath(region_overall_dir, "sensitivity_growth_$(region).png"),
        results[region].growth.figure;
        px_per_unit=DPI
    )
    save(
        joinpath(region_overall_dir, "sensitivity_survival_$(region).png"),
        results[region].survival.figure;
        px_per_unit=DPI
    )
end
