"""
Sensitivity analysis across factor space and regions
"""

include("common.jl")

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
    diameter_col::Symbol;
    n_bins::Int=10
)
    n_obs = nrow(data)
    per_bin_sample = n_obs ÷ n_bins
    bin_ids = CoralFlow.adaptive_min_sample_binning(data[!, diameter_col], per_bin_sample)

    n_bins_actual = length(unique(bin_ids))
    g_bin = Matrix{Float64}(undef, n_bins_actual, ncol(data))
    bin_details = Matrix{Float64}(undef, n_bins_actual, 3)

    for i in sort(unique(bin_ids))
        bin_sel = bin_ids .== i

        bin_mean = mean(data[bin_sel, diameter_col])
        (bin_start, bin_end) = extrema(data[bin_sel, diameter_col])
        bin_details[i, :] .= bin_start, bin_mean, bin_end

        Si = pawn(data[bin_sel, :], y_values[bin_sel]; S=10)[PAWNᵢ=At(:median)]
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

    f = Figure(; size=(800, 600))
    ax = Axis(f[1, 1])

    heatmap!(ax, g_bin)

    ax.yticks = (1:length(feature_names), feature_names)
    ax.xticks = (1:n_bins, string.(round.(bin_details[:, 2]; digits=2)))
    ax.title = "$(titlecase(replace(region, "_" => " "))) - $(analysis_type)"
    ax.ylabel = "Factors"
    ax.xlabel = "Mean Diameter of Bin\n($(per_bin_sample) samples per bin)"

    Colorbar(f[1, 2]; limits=(-0.1, max(maximum(g_bin), 1.0)), label="PAWN Index")

    return f
end

"""
    prepare_growth_data(model_results)

Prepare growth data for sensitivity analysis.
"""
function prepare_growth_data(model_results)
    all_growth = vcat(values(model_results.growth_groupings)...)
    all_y_growth = all_growth.diamnext

    ignore_cols = [g for g in growth_ignore_cols if g in propertynames(all_growth)]
    select!(all_growth, Not(ignore_cols))
    cleanup_features!(all_growth)

    return all_growth, all_y_growth
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

    ignore_cols = [g for g in surv_ignore_cols if g in propertynames(all_surv)]
    select!(all_surv, Not(ignore_cols))
    cleanup_features!(all_surv)

    return all_surv, all_y_surv
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
    model_results = CoralFlow.process_ecorrap_models(
        ecorrap_file,
        species_file;
        region=region,
        save_models=false,
        plot_validation=false,
        growth_degree=2
    )

    # Prepare data
    @info "Preparing growth data..."
    all_growth, all_y_growth = prepare_growth_data(model_results)

    @info "Preparing survival data..."
    all_surv, all_y_surv = prepare_survival_data(model_results)

    # Perform sensitivity analyses
    @info "Analyzing growth sensitivity..."
    g_bin_growth, bin_details_growth, per_bin_sample_growth = binned_sa(
        all_growth, all_y_growth, :diam; n_bins=n_bins
    )

    @info "Analyzing survival sensitivity..."
    g_bin_surv, bin_details_surv, per_bin_sample_surv = binned_sa(
        all_surv, all_y_surv, :diam_mort; n_bins=n_bins
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
regions = ["offshore_north"]  # , "torres_strait"
# EcoRRAP data for IPM_250624.csv
# ecorrap_adult_juv_combined_2021_2023_24062025
ecorrap_file = "../data/EcoRRAP data for IPM_250624.csv"
species_file = "../data/ecorrap to cscape species.csv"

results = Dict{String,NamedTuple}()

for region in regions
    results[region] = process_region_sensitivity(
        ecorrap_file, species_file, region; n_bins=10
    )
end

# Display figures
for region in regions
    @info "Displaying results for $region"
    display(results[region].growth.figure)
    display(results[region].survival.figure)
end

# Save results
for region in regions
    save(
        "$(FIG_DIR)/sensitivity/sensitivity_growth_$(region).png",
        results[region].growth.figure;
        px_per_unit=DPI
    )
    save(
        "$(FIG_DIR)/sensitivity/sensitivity_survival_$(region).png",
        results[region].survival.figure;
        px_per_unit=DPI
    )
end
