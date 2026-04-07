"""
Exploratory sensitivity analysis for the offshore north region.
"""

include("common.jl")

region_scale = ["offshore_north"]
reef_target = [nothing, "moore"]

stats_of_interest = [:mean, :median, :std]
for reg_scale in region_scale
    for rt in reef_target
        model_results = CoralFlow.process_ecorrap_models(
            "../data/EcoRRAP data for IPM_250624.csv",
            "data/ecorrap_to_cscape_species.csv";
            region=reg_scale,
            reef=rt,
            save_models=true,
            output_dir="./$(OUTPUT_DIR)/$(reg_scale)/$(isnothing(rt) ? "overall" : rt)",
            plot_validation=false,
            growth_degree=1,
            survival_degree=2,
            n_bins=10
        )

        human_region_name = titlecase(replace(reg_scale, "_" => " "))
        if !isnothing(rt)
            scale_name = reg_scale * "_" * rt
            human_reef_name = titlecase(replace(rt, "_" => " ")) * " Reef"
            human_scale_name = "$(human_reef_name) ($(human_region_name))"
        else
            scale_name = reg_scale
            human_scale_name = human_region_name
        end

        fig_output = joinpath(FIG_DIR, "sensitivity", reg_scale, isnothing(rt) ? "overall" : rt)
        mkpath(fig_output)

        ##### Overall ####

        all_growth = vcat(values(model_results.growth_groupings)...)
        all_surv = vcat(values(model_results.survival_groupings)...)

        scale_fn = replace(scale_name, " " => "_")
        reef_dir = joinpath(OUTPUT_DIR, reg_scale, isnothing(rt) ? "overall" : rt)
        mkpath(reef_dir)
        CSV.write(joinpath(reef_dir, scale_fn * "_growth.csv"), all_growth)
        CSV.write(joinpath(reef_dir, scale_fn * "_survival.csv"), all_surv)

        all_y_growth = all_growth.est_1yo_growth
        all_y_surv = all_surv.surv
        all_y_surv[ismissing.(all_y_surv), :] .= 0
        all_y_surv = Int64.(all_y_surv)

        ignore_cols = [g for g in growth_ignore_cols if g in propertynames(all_growth)]
        select!(
            all_growth,
            Not(ignore_cols)
        )

        ignore_cols = [g for g in surv_ignore_cols if g in propertynames(all_surv)]
        select!(
            all_surv,
            Not(ignore_cols)
        )

        cleanup_features!(all_growth)
        cleanup_features!(all_surv)
        rename_for_display!(all_growth)
        rename_for_display!(all_surv)

        Si_growth = pawn(all_growth, all_y_growth; S=10)[PAWNᵢ=At(stats_of_interest)]
        f, ax, sp = heatmap(
            Si_growth[sortperm(Si_growth[PAWNᵢ=At(:median)]), :];
            colorrange=(-0.1, maximum(Si_growth)), colormap=:viridis
        )
        ax.xlabelsize = 14
        ax.ylabelsize = 14
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12
        ax.titlesize = 14
        ax.xticklabelrotation = π / 8

        ax.title = "Growth - $(human_scale_name)"

        resize!(f, 800, 286)
        save("$(fig_output)/Si_$(scale_fn)_growth_overall.png", f; px_per_unit=DPI)

        # Analysis indicate that for specific locales, diameter is an influential factor.
        # But this may differ between locales, need to do further analyses.
        # At regional scales, dpeth, wave activity, size at mortality, and factors relating
        # to position matter.
        Si_surv = pawn(all_surv, convert.(Float64, all_y_surv); S=10)[PAWNᵢ=At(
            stats_of_interest
        )]
        f, ax, sp = heatmap(Si_surv[sortperm(Si_surv[PAWNᵢ=At(:median)]), :];
            colorrange=(-0.1, maximum(Si_surv)), colormap=:viridis
        )
        ax.xlabelsize = 14
        ax.ylabelsize = 14
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12
        ax.titlesize = 14
        ax.xticklabelrotation = π / 8

        ax.title = "Survival - $(human_scale_name)"

        resize!(f, 800, 286)
        save("$(fig_output)/Si_$(scale_fn)_survival_overall.png", f; px_per_unit=DPI)

        #### Group-specific analyses ####
        ####
        # But group specific analyses at specific locales indicate size/diameter is important
        # for most coral groups. While the influence of the various factors considered and their
        # rankings do differ from location to location, size is typically a common theme.
        ####

        global Si_growth_plots = []
        global Si_surv_plots = []
        global Si_growth_data = []
        global Si_surv_data = []
        for taxa in keys(model_results.growth_groupings)
            @info "Assessing $(taxa)"

            g_id = first(findall(CoralFlow.TARGET_GROUPS .== taxa))
            group_title = CoralFlow.GROUP_NAMES[g_id]

            X_growth = copy(model_results.growth_groupings[taxa])
            y_size = Float64.(model_results.growth_groupings[taxa].est_1yo_growth)

            X_surv = copy(model_results.survival_groupings[taxa])
            y_surv = try
                Int64.(model_results.survival_groupings[taxa].surv)
            catch e
                if !(e isa MethodError)  # isa(e, MethodError)
                    @info typeof(e)
                    rethrow(e)
                end

                tmp = model_results.survival_groupings[taxa]
                tmp[ismissing.(tmp.surv), :surv] .= 0
                Int64.(tmp.surv)
            end

            ignore_cols = [g for g in surv_ignore_cols if g in propertynames(X_surv)]
            select!(
                X_surv,
                Not(ignore_cols)
            )

            ignore_cols = [g for g in growth_ignore_cols if g in propertynames(X_growth)]
            select!(
                X_growth,
                Not(ignore_cols)
            )

            cleanup_features!(X_surv)
            cleanup_features!(X_growth)
            rename_for_display!(X_surv)
            rename_for_display!(X_growth)

            Si_growth = pawn(X_growth, y_size; S=10)[PAWNᵢ=At(stats_of_interest)]
            f, ax, sp = heatmap(
                Si_growth[sortperm(Si_growth[PAWNᵢ=At(:median)]), :];
                colormap=:viridis, colorrange=(-0.1, maximum(Si_growth))
            )
            ax.title = "Growth - $(human_scale_name)\n$(group_title)"
            ax.xlabelsize = 14
            ax.ylabelsize = 14
            ax.xticklabelsize = 12
            ax.yticklabelsize = 12
            ax.titlesize = 14
            ax.xticklabelrotation = π / 8

            resize!(f, 800, 286)
            push!(Si_growth_plots, f)
            push!(Si_growth_data, Si_growth)
            save("$(fig_output)/Si_$(scale_fn)_$(taxa)_growth.png", f; px_per_unit=DPI)

            Si_surv = pawn(X_surv, convert.(Float64, y_surv); S=10)[PAWNᵢ=At(
                stats_of_interest
            )]
            f, ax, sp = heatmap(
                Si_surv[sortperm(Si_surv[PAWNᵢ=At(:median)]), :];
                colormap=:viridis, colorrange=(-0.1, maximum(Si_surv))
            )
            ax.xlabelsize = 14
            ax.ylabelsize = 14
            ax.xticklabelsize = 12
            ax.yticklabelsize = 12
            ax.titlesize = 14
            ax.xticklabelrotation = π / 8

            ax.title = "Survival - $(human_scale_name)\n$(group_title)"
            resize!(f, 800, 286)
            push!(Si_surv_plots, f)
            push!(Si_surv_data, Si_surv)

            save("$(fig_output)/Si_$(scale_fn)_$(taxa)_survival.png", f; px_per_unit=DPI)
        end
    end
end

for s in Si_surv_data
    @info s.data
end
