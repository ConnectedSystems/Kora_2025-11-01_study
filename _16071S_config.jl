# Configuration for Moore Reef ensemble

# 15 x 100m EcoRRAP transect (1500m²), assume 60% is coral habitable area
# Figure 12
# https://www.aims.gov.au/sites/default/files/2023-08/AIMS_EcoRRAP_SOP14_V1_FIeld-Photogrammetry4D_Overview-infield-workflow_2023.pdf
reef_config = ReefConfig(;
    reef_id="16071S",
    reef_name="Moore Reef (16-071)",
    area=Float32(72.0 * 4),  # size of each EcoRRAP transect/plot times # of plots
    depth=9.0,  # Assumed average (12 + 5) / 2
    density=10,
    initial_proportions=[0.1f0, 0.4f0, 0.25f0, 0.05f0, 0.2f0],
    exclude_years=[2018, 2020, 2022, 2023]  # two obs after bleaching and two years of ecorrap
)

file_paths = CalibrationDataPaths(;
    dhw_scenarios="$(OUTPUT_DIR)/dhw_scens.nc",
    canonical_reefs="$(OUTPUT_DIR)/rrap_canonical_2025-07-15-T10-48-29.gpkg",
    growth_models="$(OUTPUT_DIR)/model/offshore_north_moore_growth_models.dat",
    survival_models="$(OUTPUT_DIR)/model/offshore_north_moore_survival_models.dat",
    output_dir=OUTPUT_DIR,
    figure_dir=FIG_DIR
)

# Optimization settings
opt_config = OptimizationConfig(;
    max_steps=50_000,
    population_size=50,
    fitness_threshold=0.4,
    ensemble_members=250,
    trace_interval=10,
    random_seed=78
)

# Search ranges
# Replaced by `param_bounds` further below
# search_ranges = SearchRanges(;
#     density=(2.0, 12.0),
#     group_proportion=(0.0, 1.0),  # probability levels for Gamma quantiles [0 - 1]
#     size_mean=(0.1, 5.0),
#     size_std=(0.1, 2.0),
#     scalers=(0.25, 2.0),
#     recruitment=(0.001, 0.1),
#     self_seeding=(0.001, 0.3)
# )

# Calibration settings
calib_settings = CalibrationSettings(;
    date_match_max_days=120,
    start_month=3,
    use_scalers=true
)

param_bounds = (;
    density=(2.0, 12.0),

    # Group proportions
    prop_tab_acro=(0.001, 0.6),
    prop_cor_acro=(0.001, 0.6),
    prop_cor_non_acro=(0.001, 0.6),
    prop_sm_mass=(0.001, 0.6),
    prop_lrg_mass=(0.001, 0.6),

    # Size dist mean
    size_mean_tab_acro=(0.2, 4.0),
    size_mean_cor_acro=(0.2, 4.0),
    size_mean_cor_non_acro=(0.2, 4.0),
    size_mean_sm_mass=(0.2, 4.0),
    size_mean_lrg_mass=(0.5, 5.0),

    # Size dist stdev
    size_stdev_tab_acro=(0.1, 2.5),
    size_stdev_cor_acro=(0.1, 2.5),
    size_stdev_cor_non_acro=(0.1, 2.5),
    size_stdev_sm_mass=(0.1, 2.5),
    size_stdev_lrg_mass=(0.1, 2.5),

    # Growth scalers
    scalers_tab_acro=(0.25, 1.75),
    scalers_cor_acro=(0.25, 1.75),
    scalers_cor_non_acro=(0.25, 1.75),
    scalers_sm_mass=(0.25, 1.75),
    scalers_lrg_mass=(0.25, 1.75),

    # Recruitment
    recruitment=(0.001, 0.2),
    self_seeding=(0.001, 0.3)
)
