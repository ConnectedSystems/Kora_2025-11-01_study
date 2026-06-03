include(joinpath(@__DIR__, "common.jl"))

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
    mkpath(joinpath(FIG_DIR, "regressions", region, tgt_dir))

    growth_results = CoralFlow.process_growth_models(
        joinpath(OUTPUT_DIR, "ecorrap_expanded.parquet"),
        joinpath(OUTPUT_DIR, "ecorrap_to_cscape_species.csv");
        region=region,
        reef=reef,
        output_dir=joinpath(OUTPUT_DIR, region, tgt_dir),
        degree=1,
        n_bins=10
    )

    survival_results = CoralFlow.process_survival_models(
        joinpath(OUTPUT_DIR, "ecorrap_expanded.parquet"),
        joinpath(OUTPUT_DIR, "ecorrap_to_cscape_species.csv");
        region=region,
        reef=reef,
        output_dir=joinpath(OUTPUT_DIR, region, tgt_dir),
        degree=2,
        n_bins=10
    )

    push!(region_growth, growth_results)
    push!(region_survival, survival_results)

    CoralFlow.viz.survival_performance_plots(
        survival_results.survival_groupings,
        survival_results.survival_fits;
        save_path=joinpath(FIG_DIR, "regressions", region, tgt_dir)
    )

    CoralFlow.viz.growth_performance_plots(
        growth_results.growth_groupings,
        growth_results.growth_fits;
        save_path=joinpath(FIG_DIR, "regressions", region, tgt_dir)
    )

    model_dir = joinpath(OUTPUT_DIR, region, tgt_dir)
    export_model_summaries(
        growth_results.growth_fits, model_dir, "$(region)_$(tgt_dir)_growth"
    )
    export_model_summaries(
        survival_results.survival_fits, model_dir, "$(region)_$(tgt_dir)_survival"
    )
end
