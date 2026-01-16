# Configuration for Moore Reef ensemble

# 15 x 100m EcoRRAP transect (1500m²), assume 60% is coral habitable area
# Figure 12
# https://www.aims.gov.au/sites/default/files/2023-08/AIMS_EcoRRAP_SOP14_V1_FIeld-Photogrammetry4D_Overview-infield-workflow_2023.pdf
reef_config = ReefConfig(;
    reef_id="16071S",
    reef_name="Moore Reef (16-071)",
    area=Float32(1500.0 * 0.6),
    depth=9.0,
    density=10,
    initial_proportions=[0.35f0, 0.3f0, 0.1f0, 0.1f0, 0.15f0],
    exclude_years=[2018, 2020]  # two obs after bleaching
)

file_paths = CalibrationDataPaths(;
    dhw_scenarios="$(OUTPUT_DIR)/dhw_scens.nc",
    canonical_reefs="$(OUTPUT_DIR)/rrap_canonical_2025-07-15-T10-48-29.gpkg",
    growth_models="$(OUTPUT_DIR)/offshore_north_moore_growth_models.dat",
    survival_models="$(OUTPUT_DIR)/offshore_north_moore_survival_models.dat",
    output_dir=OUTPUT_DIR,
    figure_dir=FIG_DIR
)

# Optimization settings
opt_config = OptimizationConfig(;
    max_steps=25_000,
    population_size=75,
    fitness_threshold=0.30,
    ensemble_members=1000,
    trace_interval=10,
    random_seed=78
)

# Search ranges
search_ranges = SearchRanges(;
    density=(1.0, 10.0),
    group_proportion=(0.0, 1.0),  # probability levels for Gamma quantiles
    size_mean=(1.0, 5.0),
    size_std=(0.1, 3.0),
    scalers=(0.25, 1.75),
    recruitment=(0.001, 0.2),
    self_seeding=(0.001, 0.2)
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
    prop_tab_acro=(0.001, 0.4),
    prop_cor_acro=(0.001, 0.4),
    prop_cor_non_acro=(0.001, 0.4),
    prop_sm_mass=(0.001, 0.4),
    prop_lrg_mass=(0.001, 0.4),

    # Size dist mean
    size_mean_tab_acro=(1.0, 5.0),
    size_mean_cor_acro=(1.0, 5.0),
    size_mean_cor_non_acro=(1.0, 5.0),
    size_mean_sm_mass=(1.0, 5.0),
    size_mean_lrg_mass=(1.0, 5.0),

    # Size dist stdev
    size_stdev_tab_acro=(0.1, 3.0),
    size_stdev_cor_acro=(0.1, 3.0),
    size_stdev_cor_non_acro=(0.1, 3.0),
    size_stdev_sm_mass=(0.1, 3.0),
    size_stdev_lrg_mass=(0.1, 3.0),

    # Growth scalers
    scalers_tab_acro=(0.25, 1.75),
    scalers_cor_acro=(0.25, 1.75),
    scalers_cor_non_acro=(0.25, 1.75),
    scalers_sm_mass=(0.25, 1.75),
    scalers_lrg_mass=(0.25, 1.75),

    # Recruitment
    recruitment=(0.001, 0.2),
    self_seeding=(0.001, 0.2)
)