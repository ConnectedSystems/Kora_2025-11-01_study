"""
Self-describing HDF5 storage for analysis results.

Included by `scripts/common.jl` and, standalone, by `scripts/extract_pawn_results.jl`,
so it must not depend on anything beyond HDF5 and DimensionalData.
"""

using HDF5
using DimensionalData

# ── HDF5 result I/O ───────────────────────────────────────────────────────────
#
# Analysis results are cached as self-describing HDF5 rather than Julia
# `serialize` output. `serialize` reconstructs by concrete type, so a
# dependency-level change -- as when YAXArray gave way to DimArray -- silently
# invalidates every cached file. HDF5 stores plain arrays plus enough metadata to
# rebuild the wrapper, and the files stay readable from Python and HDFView.
#
# Caching guards elsewhere key on these `.h5` paths, so a format change is also
# what forces stale results to be recomputed.
#
# Every file carries a `kind` attribute on the root group that `load_result`
# dispatches on:
#   "dimarray"   -- DimArray; dimension names and lookup values stored alongside
#   "array"      -- any array of a bits type
#   "vecvec"     -- Vector of equal-length numeric Vectors, stored as a matrix
#   "intdict"    -- Dict{Int,<:AbstractArray}, one dataset per key
#   "namedtuple" -- flat NamedTuple of the above plus scalars and string vectors
#
# Objects with no array representation (ReefState, BlackBoxOptim state, fitted
# regression models) remain on `serialize`/`.dat` -- see save_calibration_results.

"""Lookup values, plus the marker needed to restore their element type."""
function _h5_lookup(vals)
    eltype(vals) <: Symbol && return (String.(vals), "symbol")
    eltype(vals) <: AbstractString && return (String.(vals), "string")

    return (collect(vals), "numeric")
end

"""Inverse of [`_h5_lookup`](@ref): restore the element type from its marker."""
_h5_restore(vals, marker::AbstractString) = marker == "symbol" ? Symbol.(vals) : vals

"""
    save_result(path, x)

Write an analysis result to `path` as HDF5, creating parent directories as
needed. Dispatches on the type of `x`; see the section comment for the layouts.
"""
function save_result(path::String, da::AbstractDimArray)
    mkpath(dirname(path))
    h5open(path, "w") do fid
        attrs(fid)["kind"] = "dimarray"
        fid["data"] = Array(parent(da))
        attrs(fid)["dim_names"] = [string(DimensionalData.name(d)) for d in dims(da)]
        for (i, d) in enumerate(dims(da))
            vals, marker = _h5_lookup(collect(DimensionalData.lookup(d)))
            fid["dim_$(i)"] = vals
            attrs(fid)["dim_$(i)_type"] = marker
        end
    end

    return path
end
function save_result(path::String, x::AbstractArray)
    mkpath(dirname(path))
    h5open(path, "w") do fid
        attrs(fid)["kind"] = "array"
        fid["data"] = Array(x)
    end

    return path
end
function save_result(path::String, x::AbstractVector{<:AbstractVector{<:Real}})
    allequal(length.(x)) || throw(
        ArgumentError("Cannot store ragged vector-of-vectors as HDF5: $(path)")
    )

    mkpath(dirname(path))
    h5open(path, "w") do fid
        attrs(fid)["kind"] = "vecvec"
        # Columns are the original vectors, so the round trip is eachcol()
        fid["data"] = reduce(hcat, x)
    end

    return path
end
function save_result(path::String, d::AbstractDict{<:Integer,<:AbstractArray})
    mkpath(dirname(path))
    h5open(path, "w") do fid
        attrs(fid)["kind"] = "intdict"
        ks = sort(collect(keys(d)))
        attrs(fid)["keys"] = ks
        for k in ks
            fid["key_$(k)"] = Array(d[k])
        end
    end

    return path
end
function save_result(path::String, nt::NamedTuple)
    mkpath(dirname(path))
    h5open(path, "w") do fid
        attrs(fid)["kind"] = "namedtuple"
        attrs(fid)["fields"] = [string(k) for k in keys(nt)]
        for (k, v) in pairs(nt)
            name = string(k)
            if v isa AbstractVector{<:AbstractVector{<:Real}}
                allequal(length.(v)) || throw(
                    ArgumentError("Field $(name) of $(path) is ragged")
                )
                fid[name] = reduce(hcat, v)
                attrs(fid)["field_$(name)"] = "vecvec"
            elseif v isa AbstractArray
                vals, marker = _h5_lookup(collect(v))
                fid[name] = vals
                attrs(fid)["field_$(name)"] = marker == "numeric" ? "array" : marker
            else
                fid[name] = v
                attrs(fid)["field_$(name)"] = "scalar"
            end
        end
    end

    return path
end

"""
    load_result(path)

Read back a file written by [`save_result`](@ref), restoring the original type.
"""
function load_result(path::String)
    return h5open(path, "r") do fid
        kind = attrs(fid)["kind"]
        if kind == "dimarray"
            names_ = attrs(fid)["dim_names"]
            ds = ntuple(length(names_)) do i
                vals = _h5_restore(read(fid["dim_$(i)"]), attrs(fid)["dim_$(i)_type"])
                Dim{Symbol(names_[i])}(vals)
            end
            return DimArray(read(fid["data"]), ds)
        elseif kind == "array"
            return read(fid["data"])
        elseif kind == "vecvec"
            return collect.(eachcol(read(fid["data"])))
        elseif kind == "intdict"
            ks = attrs(fid)["keys"]
            return Dict(Int(k) => read(fid["key_$(k)"]) for k in ks)
        elseif kind == "namedtuple"
            fields = attrs(fid)["fields"]
            vals = map(fields) do name
                marker = attrs(fid)["field_$(name)"]
                raw = read(fid[name])
                marker == "vecvec" && return collect.(eachcol(raw))
                marker == "symbol" && return Symbol.(raw)
                return raw
            end
            return NamedTuple{Tuple(Symbol.(fields))}(Tuple(vals))
        end

        throw(ArgumentError("Unrecognised result kind '$(kind)' in $(path)"))
    end
end
