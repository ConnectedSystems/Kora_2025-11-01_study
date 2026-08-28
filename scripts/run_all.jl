"""
Run the full Kora calibration study pipeline end-to-end, in the order documented
in README.md. Run from within `scripts/`:

    julia --project=.. run_all.jl

or from the REPL (from the study root):

    ; cd scripts
    include("run_all.jl")

Each stage is wrapped so a failure is logged and the run continues to the next
stage rather than aborting the whole pipeline silently. Stages 3/4 (ensemble
calibration) are the long pole — each spawns 20 worker processes and can take
hours.
"""

const STAGES = [
    "00a_prep_ecorrap_data.jl",
    "00b_study_area_map.jl",
    "01a_exploratory_SA.jl",
    "01b_exploratory_SA_offshore_north.jl",
    "01c_exploratory_SA_torres_strait.jl",
    "01d_compile_SA_figures.jl",
    "01e_paper_dataset_summary.jl",
    "02a_offshore_north_fit_to_diameter.jl",
    "02b_torres_strait_fit_to_diameter.jl",
    "03a_moore_ensemble.jl",
    "03b_moore_ensemble_assessment.jl",
    "04a_masig_torres_strait_ensemble.jl",
    "04b_masig_ensemble_assessment.jl",
    "05_combined_sa_heatmaps.jl"
]

_PIPELINE_RESULTS = Dict{String,Symbol}()

for _pipeline_stage in STAGES
    _pipeline_path = joinpath(@__DIR__, _pipeline_stage)
    _pipeline_t0 = time()
    @info "==== Starting $_pipeline_stage ===="
    try
        include(_pipeline_path)
        _PIPELINE_RESULTS[_pipeline_stage] = :ok
        @info "==== Finished $_pipeline_stage ($(round(time() - _pipeline_t0, digits=1))s) ===="
    catch err
        _PIPELINE_RESULTS[_pipeline_stage] = :failed
        @error "==== FAILED $_pipeline_stage ($(round(time() - _pipeline_t0, digits=1))s) ====" exception=(err, catch_backtrace())
    end
end

@info "==== Pipeline summary ====" _PIPELINE_RESULTS
_pipeline_failed = [k for (k, v) in _PIPELINE_RESULTS if v == :failed]
if !isempty(_pipeline_failed)
    @warn "Stages that failed" _pipeline_failed
end
