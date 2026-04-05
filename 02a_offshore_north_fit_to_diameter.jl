include("common.jl")

region = "offshore_north"
reefs = [nothing, "moore"]
region_growth = []
region_survival = []

for reef in reefs
    if isnothing(reef)
        tgt_dir = "overall"
    else
        tgt_dir = reef
    end
    mkpath("./figs/regressions/$(region)/$(tgt_dir)")

    growth_results = CoralFlow.process_growth_models(
        "../data/ecorrap_adult_juv_combined_2021_2023_24062025.csv",
        "data/ecorrap_to_cscape_species.csv";
        region=region,
        reef=reef,
        output_dir="./$(OUTPUT_DIR)/model",
        degree=1,
        n_bins=10
    )

    survival_results = CoralFlow.process_survival_models(
        "../data/ecorrap_adult_juv_combined_2021_2023_24062025.csv",
        "data/ecorrap_to_cscape_species.csv";
        region=region,
        reef=reef,
        output_dir="./$(OUTPUT_DIR)/model",
        degree=2,
        n_bins=10
    )

    push!(region_growth, growth_results)
    push!(region_survival, survival_results)

    CoralFlow.viz.survival_performance_plots(
        survival_results.survival_groupings,
        survival_results.survival_fits;
        save_path="./figs/regressions/$(region)/$(tgt_dir)"
    )

    CoralFlow.viz.growth_performance_plots(
        growth_results.growth_groupings,
        growth_results.growth_fits;
        save_path="./figs/regressions/$(region)/$(tgt_dir)"
    )
end
