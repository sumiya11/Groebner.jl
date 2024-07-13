# This file is a part of Groebner.jl. License is GNU GPL v2.

# CompositeInt is a tuple of several integers

struct CompositeInt{N, T}
    data::NTuple{N, T}
end

function CompositeInt{N, T}(a::CompositeInt{N, U}) where {N, T, U}
    CompositeInt{N, T}(a.data .% T)
end

composite_int_unpack(a::CompositeInt{N, T}) where {N, T} = a.data

Base.typemax(::Type{CompositeInt{N, T}}) where {N, T} =
    CompositeInt(ntuple(_ -> typemax(T), N))
Base.typemin(::Type{CompositeInt{N, T}}) where {N, T} =
    CompositeInt(ntuple(_ -> typemin(T), N))

Base.rem(a::CompositeInt{N, T}, ::Type{CompositeInt{N, U}}) where {N, T, U} =
    CompositeInt(a.data .% U)

Base.widen(::Type{CompositeInt{N, T}}) where {N, T} = CompositeInt{N, widen(T)}

Base.isless(ci::CompositeInt{N, T}, cj::CompositeInt{N, U}) where {N, T, U} =
    all(ci.data .< cj.data)
Base.isless(ci::CompositeInt{N, T}, x::U) where {N, T, U <: Integer} = all(<(x), ci.data)
Base.isless(x::U, ci::CompositeInt{N, T}) where {N, T, U <: Integer} = all(x .< ci.data)

Base.iszero(ci::CompositeInt) = all(iszero, ci.data)
Base.isone(ci::CompositeInt)  = all(isone, ci.data)

Base.zero(::Type{CompositeInt{N, T}}) where {N, T} = CompositeInt(ntuple(_ -> zero(T), N))
Base.one(::Type{CompositeInt{N, T}}) where {N, T} = CompositeInt(ntuple(_ -> one(T), N))

Base.:+(a::CompositeInt{N, T}, b::CompositeInt{N, T}) where {N, T} =
    CompositeInt(a.data .+ b.data)

Base.:-(a::CompositeInt{N, T}, b::CompositeInt{N, T}) where {N, T} =
    CompositeInt(a.data .- b.data)

Base.:*(a::CompositeInt{N, T}, b::CompositeInt{N, T}) where {N, T} =
    CompositeInt(a.data .* b.data)

function Base.muladd(
    c::CompositeInt{N, T},
    a::CompositeInt{N, T},
    b::CompositeInt{N, T}
) where {N, T}
    CompositeInt(c.data .* a.data .+ b.data)
end
