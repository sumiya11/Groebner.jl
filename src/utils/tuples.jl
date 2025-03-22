# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Struct <--> NamedTuple.
# Adapted from NamedTupleTools.jl. License is MIT.

function fieldvalues(x::T) where {T}
    !isstructtype(T) && throw(ArgumentError("$(T) is not a struct type"))
    return ((getfield(x, name) for name in fieldnames(T))...,)
end

function nt_from_struct(x::T) where {T}
    !isstructtype(T) && throw(ArgumentError("$(T) is not a struct type"))
    names = fieldnames(T)
    types = fieldtypes(T)
    values = fieldvalues(x)
    return NamedTuple{names, Tuple{types...}}(values)
end

function struct_from_nt(::Type{S}, x::NT) where {S, N, T, NT <: NamedTuple{N, T}}
    names = N
    values = fieldvalues(x)
    fieldnames(S) != names && throw(ErrorException("fields in ($S) do not match ($x)"))
    return S(values...)
end

function struct_update(::Type{T}, str, new_options::NamedTuple) where {T}
    options = nt_from_struct(str)
    options = merge(options, new_options)
    struct_from_nt(T, options)
end
