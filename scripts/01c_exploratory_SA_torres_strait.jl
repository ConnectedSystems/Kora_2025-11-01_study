"""
Exploratory sensitivity analysis for the Torres Strait region.
"""

include(joinpath(@__DIR__, "common.jl"))

Random.seed!(42)

region_scale = ["torres_strait"]
reef_target = [nothing, "masig"]

stats_of_interest = [:mean, :std]
for reg_scale in region_scale
    for rt in reef_target
        model_results = Kora.process_ecorrap_models(
            joinpath(OUTPUT_DIR, "ecorrap_unified.parquet"),
            joinpath(OUTPUT_DIR, "ecorrap_to_groups.csv");
            region=reg_scale,
            reef=rt,
            save_models=true,
            output_dir=joinpath(OUTPUT_DIR, reg_scale, isnothing(rt) ? "overall" : rt),
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

        fig_output = joinpath(
            FIG_DIR, "sensitivity", reg_scale, isnothing(rt) ? "overall" : rt
        )
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

        for col in names(all_growth)
            eltype(all_growth[!, col]) <: AbstractFloat && replace!(all_growth[!, col], NaN => -1.0)
        end
        for col in names(all_surv)
            eltype(all_surv[!, col]) <: AbstractFloat && replace!(all_surv[!, col], NaN => -1.0)
        end

        Si_growth = pawn(all_growth, all_y_growth; S=10)
        f = plot_pawn_heatmap(Si_growth, "Growth - $(human_scale_name)")
        save("$(fig_output)/Si_$(scale_fn)_growth_overall.png", f; px_per_unit=DPI)

        # Analysis indicate that for specific locales, diameter is an influential factor.
        # But this may differ between locales, need to do further analyses.
        # At regional scales, depth, wave activity, size at mortality, and factors relating to position
        # matter.
        Si_surv = pawn(all_surv, convert.(Float64, all_y_surv); S=10)
        f = plot_pawn_heatmap(Si_surv, "Survival - $(human_scale_name)")
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

            g_id = first(findall(Kora.TARGET_GROUPS .== taxa))
            group_title = Kora.GROUP_NAMES[g_id]

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

            for col in names(X_growth)
                eltype(X_growth[!, col]) <: AbstractFloat && replace!(X_growth[!, col], NaN => -1.0)
            end
            for col in names(X_surv)
                eltype(X_surv[!, col]) <: AbstractFloat && replace!(X_surv[!, col], NaN => -1.0)
            end

            Si_growth = pawn(X_growth, y_size; S=10)
            f = plot_pawn_heatmap(Si_growth, "Growth - $(human_scale_name)\n$(group_title)")
            push!(Si_growth_plots, f)
            push!(Si_growth_data, Si_growth)
            save("$(fig_output)/Si_$(scale_fn)_$(taxa)_growth.png", f; px_per_unit=DPI)

            Si_surv = pawn(X_surv, convert.(Float64, y_surv); S=10)
            f = plot_pawn_heatmap(Si_surv, "Survival - $(human_scale_name)\n$(group_title)")
            push!(Si_surv_plots, f)
            push!(Si_surv_data, Si_surv)
            save("$(fig_output)/Si_$(scale_fn)_$(taxa)_survival.png", f; px_per_unit=DPI)
        end
    end
end

for s in Si_surv_data
    @info s.data
end

# ─── Wave-subset sensitivity analysis ────────────────────────────────────────
#
# wave_hs_mean is sparsely covered (good for TSMA D; ~1 year for ONMO D) so it
# is excluded from the main SA to preserve full sample size.  Here we repeat
# the overall SA on the subset of observations that DO have wave data, with
# wave_hs_mean included as a feature.  The reduced-n runs are saved separately
# so figures from the main analysis are not overwritten.
#
# Comparison between the wave-subset run (with wave) and the same subset run
# (without wave) reveals wave's marginal contribution independent of
# sample-size changes.

@info "Running wave-subset SA for Torres Strait"

for reg_scale in region_scale
    for rt in reef_target
        # Reload the collated groupings written above
        scale_fn   = isnothing(rt) ? reg_scale : reg_scale * "_" * rt
        reef_dir   = joinpath(OUTPUT_DIR, reg_scale, isnothing(rt) ? "overall" : rt)
        wave_fig_dir = joinpath(
            FIG_DIR, "sensitivity", reg_scale,
            isnothing(rt) ? "overall" : rt, "wave_subset"
        )
        mkpath(wave_fig_dir)

        human_region_name = titlecase(replace(reg_scale, "_" => " "))
        human_scale_name = if isnothing(rt)
            human_region_name
        else
            titlecase(replace(rt, "_" => " ")) * " Reef ($(human_region_name))"
        end

        all_growth = CSV.read(joinpath(reef_dir, scale_fn * "_growth.csv"), DataFrame)
        all_surv   = CSV.read(joinpath(reef_dir, scale_fn * "_survival.csv"), DataFrame)

        # Filter to rows with wave data
        wave_growth = filter(r -> !ismissing(r.wave_hs_mean) && !isnan(r.wave_hs_mean), all_growth)
        wave_surv   = filter(r -> !ismissing(r.wave_hs_mean) && !isnan(r.wave_hs_mean), all_surv)

        @info "Wave subset sizes" growth=nrow(wave_growth) survival=nrow(wave_surv) scale=scale_fn

        if nrow(wave_growth) == 0 || nrow(wave_surv) == 0
            @warn "No wave data for $(scale_fn), skipping wave-subset SA"
            continue
        end

        # ── helper: build SA-ready feature matrix with wave optionally included ──
        function prepare_wave_features(df::DataFrame, ignore_cols_base, include_wave::Bool)
            # Always drop the redundant wave columns; keep wave_hs_mean only when requested.
            wave_always_drop = [:wave_hs_median, :n_days_waves]
            wave_keep = include_wave ? [:wave_hs_mean] : Symbol[]
            # Build the drop list from the base ignore list, stripping any columns we
            # want to keep, then add the always-drop wave columns.
            base_drop = filter(c -> c ∉ wave_keep, ignore_cols_base)
            drop = [c for c in vcat(base_drop, wave_always_drop) if c in propertynames(df)]
            X = select(df, Not(unique(drop)))
            cleanup_features!(X)
            rename_for_display!(X)
            return X
        end

        for include_wave in (false, true)
            suffix = include_wave ? "with_wave" : "no_wave"
            title_suffix = include_wave ? " [wave incl.]" : " [wave excl.]"

            # ── growth ──
            y_growth = wave_growth.est_1yo_growth
            X_growth = prepare_wave_features(wave_growth, growth_ignore_cols, include_wave)
            for col in names(X_growth)
                eltype(X_growth[!, col]) <: AbstractFloat && replace!(X_growth[!, col], NaN => -1.0)
            end

            Si_g = pawn(X_growth, y_growth; S=10)
            f = plot_pawn_heatmap(
                Si_g,
                "Growth (wave subset) - $(human_scale_name)$(title_suffix)\nn = $(nrow(X_growth))"
            )
            save(joinpath(wave_fig_dir, "Si_$(scale_fn)_growth_$(suffix).png"), f; px_per_unit=DPI)

            # ── survival ──
            y_surv_raw = wave_surv.surv
            y_surv = Int64.(coalesce.(y_surv_raw, 0))
            X_surv = prepare_wave_features(wave_surv, surv_ignore_cols, include_wave)
            for col in names(X_surv)
                eltype(X_surv[!, col]) <: AbstractFloat && replace!(X_surv[!, col], NaN => -1.0)
            end

            Si_s = pawn(X_surv, convert.(Float64, y_surv); S=10)
            f = plot_pawn_heatmap(
                Si_s,
                "Survival (wave subset) - $(human_scale_name)$(title_suffix)\nn = $(nrow(X_surv))"
            )
            save(joinpath(wave_fig_dir, "Si_$(scale_fn)_survival_$(suffix).png"), f; px_per_unit=DPI)
        end
    end
end
