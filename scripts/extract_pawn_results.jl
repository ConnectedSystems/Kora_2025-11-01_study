"""
Extract PAWN sensitivity results from serialized .dat files and save as CSV.

Run from within scripts/:
    julia --project=.. extract_pawn_results.jl
"""

using Serialization
using CSV, DataFrames
using DimensionalData

# NOTE: `*_pawn_results.dat` written before the DimensionalData migration hold `YAXArray`s
# and can no longer be deserialized. Re-run the corresponding `*_ensemble_assessment.jl`
# script to regenerate them as `DimArray`s.

OUTPUT_DIR = joinpath(@__DIR__, "..", "data")

"""
    pawn_to_dataframe(results::AbstractDimArray) -> DataFrame

Convert a PAWN result cube (factors × PAWNᵢ stats) to a tidy DataFrame.
"""
function pawn_to_dataframe(results::AbstractDimArray)
    factors = string.(collect(dims(results, 1)))
    stats   = collect(dims(results, 2))
    df = DataFrame(:parameter => factors)
    for s in stats
        df[!, string(s)] = collect(results[PAWNᵢ=At(s)])
    end
    return df
end

configs = [
    (
        label    = "moore (16071S)",
        reef_id  = "16071S",
        data_dir = joinpath(OUTPUT_DIR, "sensitivity", "offshore_north", "moore", "ensemble"),
    ),
    (
        label    = "masig",
        reef_id  = "masig",
        data_dir = joinpath(OUTPUT_DIR, "sensitivity", "torres_strait", "masig", "ensemble"),
    ),
]

for cfg in configs
    @info "Processing $(cfg.label)"

    for kind in ("unconstrained", "constrained")
        fn = joinpath(cfg.data_dir, "$(cfg.reef_id)_$(kind)_pawn_results.dat")
        if !isfile(fn)
            @warn "  Missing: $fn — skipping"
            continue
        end

        results = deserialize(fn)
        df = pawn_to_dataframe(results)

        out_path = joinpath(cfg.data_dir, "$(cfg.reef_id)_$(kind)_pawn_results.csv")
        CSV.write(out_path, df)
        @info "  Saved $kind → $out_path"

        # Print a ranked summary (median, descending)
        sorted = sort(df, :mean; rev=true)
        println("\n  $(cfg.label) — $(kind) PAWN (ranked by mean):")
        println("  ", rpad("Parameter", 28), rpad("mean", 9), rpad("median", 9), "cv")
        for row in eachrow(sorted)
            println(
                "  ",
                rpad(row.parameter, 28),
                rpad(round(row.mean;   digits=3), 9),
                rpad(round(row.median; digits=3), 9),
                round(row.cv; digits=3)
            )
        end
        println()
    end
end
