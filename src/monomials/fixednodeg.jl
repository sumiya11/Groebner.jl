using SmallCollections: FixedVector, SmallVector, fixedvector, sum_fast, bits, bitsize

struct FixedMonomNoDeg{N,T}
    ev::FixedVector{N,T}
end

monom_max_vars(a::FixedMonomNoDeg) = monom_max_vars(typeof(a))

monom_max_vars(::Type{<:FixedMonomNoDeg{N}}) where N = N
monom_totaldeg(a::FixedMonomNoDeg) = reduce(+, a.ev; init = zero(UInt16)) % UInt
monom_copy(a::FixedMonomNoDeg) = a
monom_entrytype(::M) where M <: FixedMonomNoDeg = monom_entrytype(M)
monom_entrytype(::Type{FixedMonomNoDeg{N,T}}) where {N,T} = T

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{<:FixedMonomNoDeg{N}}, ::Integer) where N
    rand(rng, FixedVector{N,MonomHash})
end

function monom_construct_from_vector(::Type{FixedMonomNoDeg{N,T}}, ev::AbstractVector) where {N,T}
    all(<=(typemax(T)), ev) || __throw_monom_overflow_error()
    FixedMonomNoDeg(fixedvector(SmallVector{N,T}(ev)))
end

function monom_construct_const(::Type{FixedMonomNoDeg{N,T}}, ::Integer) where {N,T}
    FixedMonomNoDeg(zero(FixedVector{N,T}))
end

function monom_hash(a::FixedMonomNoDeg{N}, b::FixedVector{N}) where N
    sum_fast(map(*, a.ev, b)) % MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::FixedMonomNoDeg)
    copyto!(tmp, 1, a.ev, 1, length(tmp))
end

function monom_is_supported_ordering(::Type{<:FixedMonomNoDeg}, ::Ord) where Ord
    Ord <: Union{DegLex{true}, DegRevLex{true}, InputOrdering}
end

function monom_isless(a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}, ::DegLex{true}) where N
    adeg = monom_totaldeg(a)
    bdeg = monom_totaldeg(b)
    adeg < bdeg || (adeg == bdeg && a.ev < b.ev)
end

function monom_isless(a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}, ::DegRevLex{true}) where N
    adeg = monom_totaldeg(a)
    bdeg = monom_totaldeg(b)
    if bitsize(a) <= 512
        adeg < bdeg || (adeg == bdeg && bits(a.ev) > bits(b.ev))
    else
        adeg != bdeg && return adeg < bdeg
        aev = reinterpret(UInt64, a.ev)
        bev = reinterpret(UInt64, b.ev)
        @inbounds for i in length(aev):-1:1
            aev[i] != bev[i] && return aev[i] > bev[i]
        end
        false
    end
end

function monom_lcm!(_, a::FixedMonomNoDeg{N,T}, b::FixedMonomNoDeg{N,T}) where {N,T}
    FixedMonomNoDeg(map(max, a.ev, b.ev))
end

function monom_lcm!(_, a::FixedMonomNoDeg{8,UInt8}, b::FixedMonomNoDeg{8,UInt8})
    # conversion to UInt16 to force vectorization
    ev = max.(a.ev .% UInt16, b.ev .% UInt16)
    FixedMonomNoDeg(ev .% UInt8)
end

function monom_is_gcd_const(a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}) where N
    iszero(map(min, a.ev, b.ev))
end

function monom_product!(_, a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}) where N
    ev, s = Base.Checked.add_with_overflow(a.ev, b.ev)
    isempty(s) || __throw_monom_overflow_error()
    FixedMonomNoDeg(ev)
end

function monom_division!(_, a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}) where N
    FixedMonomNoDeg(a.ev - b.ev)
end

function monom_is_divisible(a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}) where N
    all(a.ev .>= b.ev)
end

function monom_is_divisible!(c::FixedMonomNoDeg{N}, a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}) where N
    ev, s = Base.Checked.sub_with_overflow(a.ev, b.ev)
    isempty(s), FixedMonomNoDeg(ev)
end

function monom_is_equal(a::FixedMonomNoDeg{N}, b::FixedMonomNoDeg{N}) where N
    a.ev == b.ev
end

function monom_create_divmask(
    a::FixedMonomNoDeg{N,T},
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
