# This file is a part of Groebner.jl. License is GNU GPL v2.

# CompositeNumber is a tuple of several numbers

struct CompositeNumber{N, T}
    data::NTuple{N, T}
end

function CompositeNumber{N, T}(a::CompositeNumber{N, U}) where {N, T, U}
    CompositeNumber{N, T}(a.data .% T)
end

Base.isinteger(ci::CompositeNumber{N, T}) where {N, T} = all(isinteger.(ci.data))

Base.typemax(::Type{CompositeNumber{N, T}}) where {N, T} =
    CompositeNumber(ntuple(_ -> typemax(T), N))
Base.typemin(::Type{CompositeNumber{N, T}}) where {N, T} =
    CompositeNumber(ntuple(_ -> typemin(T), N))

Base.rem(a::CompositeNumber{N, T}, ::Type{CompositeNumber{N, U}}) where {N, T, U} =
    CompositeNumber(a.data .% U)

Base.widen(::Type{CompositeNumber{N, T}}) where {N, T} = CompositeNumber{N, widen(T)}

Base.isless(ci::CompositeNumber{N, T}, cj::CompositeNumber{N, U}) where {N, T, U} =
    all(ci.data .< cj.data)
Base.isless(ci::CompositeNumber{N, T}, x::U) where {N, T, U <: Number} = all(<(x), ci.data)
Base.isless(x::U, ci::CompositeNumber{N, T}) where {N, T, U <: Number} = all(x .< ci.data)

Base.iszero(ci::CompositeNumber) = all(iszero, ci.data)
Base.isone(ci::CompositeNumber)  = all(isone, ci.data)

Base.zero(::Type{CompositeNumber{N, T}}) where {N, T} =
    CompositeNumber(ntuple(_ -> zero(T), N))
Base.one(::Type{CompositeNumber{N, T}}) where {N, T} =
    CompositeNumber(ntuple(_ -> one(T), N))

Base.:+(a::CompositeNumber{N, T}, b::CompositeNumber{N, T}) where {N, T} =
    CompositeNumber(a.data .+ b.data)

Base.:-(a::CompositeNumber{N, T}, b::CompositeNumber{N, T}) where {N, T} =
    CompositeNumber(a.data .- b.data)

Base.:*(a::CompositeNumber{N, T}, b::CompositeNumber{N, T}) where {N, T} =
    CompositeNumber(a.data .* b.data)

function Base.muladd(
    c::CompositeNumber{N, T},
    a::CompositeNumber{N, T},
    b::CompositeNumber{N, T}
) where {N, T}
    CompositeNumber(c.data .* a.data .+ b.data)
end
