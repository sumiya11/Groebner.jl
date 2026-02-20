using SmallCollections: FixedVector, SmallVector, fixedvector, sum_fast, bits

struct FixedMonom{N,T}
    ev::FixedVector{N,T}
    deg::T
end

function FixedMonom(ev::FixedVector)
    a = FixedMonom(ev, sum_fast(ev))
    monom_overflow_check(a)
    a
end

monom_max_vars(a::FixedMonom) = monom_max_vars(typeof(a))

function monom_overflow_check(a::FixedMonom{N,T}) where {N,T}
    signbit(signed(a.deg)) && __throw_monom_overflow_error(a.deg, T)
end

monom_max_vars(::Type{<:FixedMonom{N}}) where N = N
monom_totaldeg(a::FixedMonom) = a.deg % UInt
monom_copy(a::FixedMonom) = a
monom_entrytype(::M) where M <: FixedMonom = monom_entrytype(M)
monom_entrytype(::Type{FixedMonom{N,T}}) where {N,T} = T

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{<:FixedMonom{N}}, ::Integer) where N
    rand(rng, FixedVector{N,MonomHash})
end

function monom_construct_from_vector(::Type{FixedMonom{N,T}}, ev::AbstractVector) where {N,T}
    FixedMonom(fixedvector(SmallVector{N,T}(ev)))
end

function monom_construct_const(::Type{FixedMonom{N,T}}, ::Integer) where {N,T}
    FixedMonom(zero(FixedVector{N,T}))
end

function monom_hash(a::FixedMonom{N}, b::FixedVector{N}) where N
    sum_fast(map(*, a.ev, b)) % MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::FixedMonom)
    copyto!(tmp, 1, a.ev, 1, length(tmp))
end

function monom_is_supported_ordering(::Type{<:FixedMonom}, ::Ord) where Ord
    Ord <: Union{DegLex{true}, DegRevLex{true}, InputOrdering}
end

function monom_isless(a::FixedMonom{N}, b::FixedMonom{N}, ::DegLex{true}) where N
    a.deg < b.deg || (a.deg == b.deg && a.ev < b.ev)
end

function monom_isless(a::FixedMonom{N}, b::FixedMonom{N}, ::DegRevLex{true}) where N
    a.deg < b.deg || (a.deg == b.deg && bits(a.ev) > bits(b.ev))
end

function monom_lcm!(_, a::FixedMonom{N,T}, b::FixedMonom{N,T}) where {N,T}
    FixedMonom(map(max, a.ev, b.ev))
end

function monom_is_gcd_const(a::FixedMonom{N}, b::FixedMonom{N}) where N
    iszero(bits(map(!iszero, a.ev)) & bits(map(!iszero, b.ev)))
end

function monom_product!(_, a::FixedMonom{N}, b::FixedMonom{N}) where N
    c = FixedMonom(a.ev + b.ev, a.deg + b.deg)
    monom_overflow_check(c)
    c
end

function monom_division!(_, a::FixedMonom{N}, b::FixedMonom{N}) where N
    FixedMonom(a.ev - b.ev, a.deg - b.deg)
end

function monom_is_divisible(a::FixedMonom{N}, b::FixedMonom{N}) where N
    all(map(>=, a.ev, b.ev))
end

function monom_is_divisible!(c::FixedMonom{N}, a::FixedMonom{N}, b::FixedMonom{N}) where N
    monom_is_divisible(a, b) ? (true, monom_division!(c, a, b)) : (false, c)
end

function monom_is_equal(a::FixedMonom{N}, b::FixedMonom{N}) where N
    a.ev == b.ev
end

function monom_create_divmask(
    a::FixedMonom{N,T},
    ::Type{Mask},
    ndivvars::Int,
    divmap::Vector{U},
    ndivbits::Int,
    strategy::Symbol
) where {N, T, Mask, U}
    @invariant strategy == :first_variables
    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    @inbounds for i in 1:ndivvars
        for _ in 1:ndivbits
            if a.ev[i] >= divmap[ctr]
                res |= o << ((ctr - 1) % UInt)
            end
            ctr += o
        end
    end
    res
end
