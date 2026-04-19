# Configuration for Masig Reef ensemble

reef_config = ReefConfig(;
    reef_id="masig",
    reef_name="Masig Reef",
    area=Float32(72.0 * 4),
    depth=9.0,
    density=10,
    initial_proportions=[0.02f0, 0.18f0, 0.2f0, 0.3f0, 0.3f0],
    exclude_years=Int[2024],  # last year held out for assessment
    disturbance_years=Int[]   # update if known disturbance years apply
)

file_paths = CalibrationDataPaths(;
    dhw_scenarios=joinpath(OUTPUT_DIR, "dhw_scens.nc"),
    canonical_reefs=joinpath(OUTPUT_DIR, "rrap_canonical_2025-07-15-T10-48-29.gpkg"),
    growth_models=joinpath(OUTPUT_DIR, "torres_strait", "masig", "torres_strait_masig_growth_models.dat"),
    survival_models=joinpath(OUTPUT_DIR, "torres_strait", "masig", "torres_strait_masig_survival_models.dat"),
    output_dir=OUTPUT_DIR,
    figure_dir=FIG_DIR
)

# Optimization settings
opt_config = OptimizationConfig(;
    max_steps=50_000,
    population_size=50,
    fitness_threshold=0.8,
    ensemble_members=250,
    trace_interval=10,
    random_seed=64
)

# Calibration settings
calib_settings = CalibrationSettings(;
    date_match_max_days=120,
    start_month=3,
    use_scalers=true
)

param_bounds = (;
    density=(1.0, 10.0),

    # Group proportions
    prop_tab_acro=(0.001, 0.6),
    prop_cor_acro=(0.001, 0.6),
    prop_cor_non_acro=(0.001, 0.6),
    prop_sm_mass=(0.001, 0.6),
    prop_lrg_mass=(0.001, 0.6),

    # Size dist mean
    size_mean_tab_acro=(1.0, 5.0),
    size_mean_cor_acro=(1.0, 5.0),
    size_mean_cor_non_acro=(1.0, 5.0),
    size_mean_sm_mass=(1.0, 5.0),
    size_mean_lrg_mass=(1.0, 5.0),

    # Size dist stdev
    size_stdev_tab_acro=(0.25, 2.0),
    size_stdev_cor_acro=(0.25, 2.0),
    size_stdev_cor_non_acro=(0.25, 2.0),
    size_stdev_sm_mass=(0.25, 2.0),
    size_stdev_lrg_mass=(0.25, 2.0),

    # Growth scalers — kept near 1.0 because EcoRRAP-derived growth functions
    # already encode reef-specific rates. Allowing large deviations risks finding
    # solutions that fit the calibration period under suppressive DHW but then
    # accelerate unrealistically once DHW drops.
    scalers_tab_acro=(0.95, 1.05),
    scalers_cor_acro=(0.95, 1.05),
    scalers_cor_non_acro=(0.95, 1.05),
    scalers_sm_mass=(0.95, 1.05),
    scalers_lrg_mass=(0.95, 1.05),

    # Recruitment
    recruitment=(0.001, 0.15),
    self_seeding=(0.001, 0.3)
)
