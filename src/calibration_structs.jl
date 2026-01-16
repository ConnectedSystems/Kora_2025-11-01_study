"""
    ReefConfig

Configuration for a specific reef being calibrated.

# Fields
- `reef_id::String`: Unique reef identifier (e.g., "16071S")
- `reef_name::String`: Human-readable reef name
- `area::Float32`: Reef area in m²
- `depth::Float32`: Mean depth in meters
- `density::Int`: Initial population density per m²
- `initial_proportions::Vector{Float32}`: Initial proportion of each functional group
- `exclude_years::Vector{Int}`: Years to exclude from calibration (e.g., cyclone years)
"""
struct ReefConfig
    reef_id::String
    reef_name::String
    area::Float32
    depth::Float64
    density::Int
    initial_proportions::Vector{Float32}
    exclude_years::Vector{Int}
end

"""
    ReefConfig(; kwargs...)

Create reef configuration with sensible defaults.
"""
function ReefConfig(;
    reef_id::String,
    reef_name::String=reef_id,
    area::Float32=900.0f0,
    depth::Float64=7.0,
    density::Int=10,
    initial_proportions::Vector{Float32}=[0.35f0, 0.3f0, 0.1f0, 0.1f0, 0.15f0],
    exclude_years::Vector{Int}=Int[]
)
    return ReefConfig(
        reef_id, reef_name, area, depth, density,
        initial_proportions, exclude_years
    )
end

"""
    FilePaths

Paths to input and output files for calibration.

# Fields
- `dhw_scenarios::String`: Path to DHW NetCDF file
- `canonical_reefs::String`: Path to canonical reefs GeoPackage
- `growth_models::String`: Path to serialized growth models
- `survival_models::String`: Path to serialized survival models
- `output_dir::String`: Directory for saving results
- `figure_dir::String`: Directory for saving figures
"""
struct CalibrationDataPaths
    dhw_scenarios::String
    canonical_reefs::String
    growth_models::String
    survival_models::String
    output_dir::String
    figure_dir::String
end

"""
    FilePaths(; kwargs...)

Create file paths configuration with typical defaults.
"""
function CalibrationDataPaths(;
    dhw_scenarios::String="$(OUTPUT_DIR)/dhw_scens.nc",
    canonical_reefs::String="$(OUTPUT_DIR)/rrap_canonical_2025-07-15-T10-48-29.gpkg",
    growth_models::String="$(OUTPUT_DIR)/models/offshore_north_moore_growth_models.dat",
    survival_models::String="$(OUTPUT_DIR)/models/offshore_north_moore_survival_models.dat",
    output_dir::String="data",
    figure_dir::String="figures"
)
    return CalibrationDataPaths(
        dhw_scenarios, canonical_reefs,
        growth_models, survival_models, output_dir, figure_dir
    )
end

"""
    OptimizationConfig

Configuration for the optimization algorithm.

# Fields
- `max_steps::Int`: Maximum optimization steps
- `population_size::Int`: Size of optimization population
- `fitness_threshold::Float64`: Threshold for tracking good candidates
- `trace_interval::Int`: How often to print optimization progress
- `random_seed::Int`: Random seed for reproducibility
"""
struct OptimizationConfig
    max_steps::Int
    population_size::Int
    fitness_threshold::Float64
    ensemble_members::Int
    trace_interval::Int
    random_seed::Int
end

"""
    OptimizationConfig(; kwargs...)

Create optimization configuration with default values.
"""
function OptimizationConfig(;
    max_steps::Int=50_000,
    population_size::Int=50,
    fitness_threshold::Float64=0.3,
    ensemble_members::Int=100,
    trace_interval::Int=10,
    random_seed::Int=76
)
    return OptimizationConfig(
        max_steps, population_size, fitness_threshold, ensemble_members,
        trace_interval, random_seed
    )
end

"""
    SearchRanges

Parameter search ranges for optimization.

# Fields
- `density::Tuple{Float64,Float64}`: Min/max population density per m²
- `group_proportion::Tuple{Float64,Float64}`: Min/max proportion for each group
- `size_mean::Tuple{Float64,Float64}`: Min/max for lognormal mean
- `size_std::Tuple{Float64,Float64}`: Min/max for lognormal stdev
- `scalers::Tuple{Float64,Float64}`: Min/max for growth/survival scalers
- `recruitment::Tuple{Float64,Float64}`: Min/max recruitment proportion
- `self_seeding::Tuple{Float64,Float64}`: Min/max self-seeding proportion
"""
struct SearchRanges
    density::Tuple{Float64,Float64}
    group_proportion::Tuple{Float64,Float64}
    size_mean::Tuple{Float64,Float64}
    size_std::Tuple{Float64,Float64}
    scalers::Tuple{Float64,Float64}
    recruitment::Tuple{Float64,Float64}
    self_seeding::Tuple{Float64,Float64}
end

"""
    SearchRanges(; kwargs...)

Create search ranges with default values.
"""
function SearchRanges(;
    density::Tuple{Float64,Float64}=(1.0, 10.0),
    group_proportion::Tuple{Float64,Float64}=(0.001, 0.4),
    size_mean::Tuple{Float64,Float64}=(1.0, 5.0),
    size_std::Tuple{Float64,Float64}=(0.25, 1.75),
    scalers::Tuple{Float64,Float64}=(0.5, 2.0),
    recruitment::Tuple{Float64,Float64}=(0.001, 0.1),
    self_seeding::Tuple{Float64,Float64}=(0.001, 0.3)
)
    return SearchRanges(
        density, group_proportion, size_mean, size_std,
        scalers, recruitment, self_seeding
    )
end

"""
    CalibrationSettings

General calibration settings.

# Fields
- `date_match_max_days::Int`: Maximum days difference for matching dates
- `start_month::Int`: Month to start simulation (1-12)
- `use_scalers::Bool`: Whether to optimize growth/survival scalers
"""
struct CalibrationSettings
    date_match_max_days::Int
    start_month::Int
    use_scalers::Bool
end

"""
    CalibrationSettings(; kwargs...)

Create calibration settings with default values.
"""
function CalibrationSettings(;
    date_match_max_days::Int=120,
    start_month::Int=3,
    use_scalers::Bool=true
)
    return CalibrationSettings(date_match_max_days, start_month, use_scalers)
end
