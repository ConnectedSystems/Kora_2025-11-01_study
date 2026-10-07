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
    growth_models=joinpath(OUTPUT_DIR, "torres_strait", "masig", "torres_strait_masig_growth_models.json"),
    survival_models=joinpath(OUTPUT_DIR, "torres_strait", "masig", "torres_strait_masig_survival_models.json"),
    output_dir=OUTPUT_DIR,
    figure_dir=FIG_DIR
)

# Optimization settings
opt_config = OptimizationConfig(;
    max_steps=150_000,
    population_size=100,
    fitness_threshold=0.3,
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

    # Group proportions — narrowed from a uniform (0.001, 0.6) for every group to
    # per-group ranges centred on the observed 2021-2024 ecorrap_benthic composition
    # (data/ecorrap_benthic/masig_estimate.csv), which is stable across all 4 years:
    # tab_acro ~2-3%, cor_acro ~9-10%, cor_non_acro ~6-8%, sm_mass ~10-15%,
    # lrg_mass ~67-69% of cover. With only 3 calibration points and 23 free
    # parameters, the old uniform bounds (which didn't even permit lrg_mass's true
    # ~0.68 share, capped at 0.6) let the optimizer trade off unrealistic group
    # compositions against each other to fit calibration noise.
    prop_tab_acro=(0.001, 0.08),
    prop_cor_acro=(0.02, 0.20),
    prop_cor_non_acro=(0.01, 0.15),
    prop_sm_mass=(0.05, 0.25),
    prop_lrg_mass=(0.40, 0.85),

    # Size dist mean — μ of the truncated log-normal initial colony diameter (cm,
    # see `set_population!`/`initialize_coral_population!`). Narrowed from a
    # uniform (1.0, 5.0) for every group to per-group ranges derived from
    # log(diam_t1_cm) in data/torres_strait/masig/torres_strait_masig_growth_fitdata.csv
    # (mean ± ~1.5 SD): tab_acro 1.4±1.03, cor_acro 2.56±0.99, cor_non_acro 2.72±0.74,
    # sm_mass 1.68±0.98, lrg_mass 2.67±0.83.
    size_mean_tab_acro=(0.3, 2.9),
    size_mean_cor_acro=(1.1, 4.0),
    size_mean_cor_non_acro=(1.6, 3.8),
    size_mean_sm_mass=(0.3, 3.1),
    size_mean_lrg_mass=(1.4, 3.9),

    # Size dist stdev — narrowed from a uniform (0.25, 2.0) to (0.4, 1.4), bracketing
    # the observed log(diam) SD across all 5 groups (0.74-1.03) in the same fit data.
    size_stdev_tab_acro=(0.4, 1.4),
    size_stdev_cor_acro=(0.4, 1.4),
    size_stdev_cor_non_acro=(0.4, 1.4),
    size_stdev_sm_mass=(0.4, 1.4),
    size_stdev_lrg_mass=(0.4, 1.4),

    # Growth scalers — fixed near 1.0 (not calibrated). EcoRRAP-derived growth
    # functions already encode reef-specific rates, so allowing deviation risks
    # solutions that fit the calibration period under suppressive DHW but then
    # accelerate unrealistically once DHW drops. Bounds can't be an exact point
    # (1.0, 1.0): uniform sampling and the nearest-neighbour diversity
    # normalization both require nonzero width.
    scalers_tab_acro=(0.999, 1.001),
    scalers_cor_acro=(0.999, 1.001),
    scalers_cor_non_acro=(0.999, 1.001),
    scalers_sm_mass=(0.999, 1.001),
    scalers_lrg_mass=(0.999, 1.001),

    # Recruitment
    recruitment=(0.001, 0.15),
    self_seeding=(0.001, 0.3)
)
