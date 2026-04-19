using SmallCollections: FixedVector, SmallVector, fixedvector, sum_fast, dot_fast, bits

struct FixedMonom{N, T}
    ev::FixedVector{N, T}
    deg::T
end

function FixedMonom(ev::FixedVector{N, T}) where {N, T}
    @invariant sum(Int, ev) <= typemax(UInt16)
    d = sum_fast(ev; init=zero(UInt16))
    d > typemax(T) && __throw_monom_overflow_error()
    FixedMonom(ev, d % T)
end

monom_max_vars(a::FixedMonom) = monom_max_vars(typeof(a))

monom_max_vars(::Type{<:FixedMonom{N}}) where {N} = N
monom_totaldeg(a::FixedMonom{N, T}) where {N, T} = a.deg % T
monom_copy(a::FixedMonom) = a

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{<:FixedMonom{N}}, ::Integer) where {N}
    rand(rng, FixedVector{N, MonomHash})
end

function monom_construct_from_vector(::Type{FixedMonom{N, T}}, ev::AbstractVector) where {N, T}
    all(<=(typemax(T)), ev) || __throw_monom_overflow_error()
    FixedMonom(fixedvector(SmallVector{N, T}(ev)))
end

function monom_construct_const(::Type{FixedMonom{N, T}}, ::Integer) where {N, T}
    FixedMonom(zero(FixedVector{N, T}))
end

function monom_hash(a::FixedMonom{N}, b::FixedVector{N}) where {N}
    dot_fast(a.ev, b)::MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::FixedMonom)
    copyto!(tmp, 1, a.ev, 1, length(tmp))
end

function monom_is_supported_ordering(::Type{<:FixedMonom}, ::Ord) where {Ord}
    Ord <: Union{DegLex{true}, DegRevLex{true}, InputOrdering}
end

function monom_isless(a::FixedMonom{N}, b::FixedMonom{N}, ::DegLex{true}) where {N}
    a.deg < b.deg || (a.deg == b.deg && a.ev < b.ev)
end

function monom_isless(a::FixedMonom{N}, b::FixedMonom{N}, ::DegRevLex{true}) where {N}
    if sizeof(a) <= 64
        a.deg < b.deg || (a.deg == b.deg && bits(a.ev) > bits(b.ev))
    else
        a.deg != b.deg && return a.deg < b.deg
        aev = reinterpret(UInt64, a.ev)
        bev = reinterpret(UInt64, b.ev)
        @inbounds for i in length(aev):-1:1
            aev[i] != bev[i] && return aev[i] > bev[i]
        end
        false
    end
end

function monom_lcm!(_, a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    ev = max.(a.ev, b.ev)
    d = sum_fast(ev)
    d < a.deg && __throw_monom_overflow_error()
    FixedMonom(ev, d)
end

function monom_lcm!(_, a::FixedMonom{8, UInt8}, b::FixedMonom{8, UInt8})
    # conversion to UInt16 to force vectorization
    ev = max.(a.ev .% UInt16, b.ev .% UInt16)
    d = sum_fast(ev) % UInt8
    d < a.deg && __throw_monom_overflow_error()
    FixedMonom(ev .% UInt8, d)
end

function monom_is_gcd_const(a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    iszero(min.(a.ev, b.ev))
end

function monom_product!(_, a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    d = a.deg + b.deg
    d < a.deg && __throw_monom_overflow_error()
    FixedMonom(a.ev + b.ev, d)
end

function monom_division!(_, a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    FixedMonom(a.ev - b.ev, a.deg - b.deg)
end

function monom_is_divisible(a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    all(a.ev .>= b.ev)
end

function monom_is_divisible!(_, a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    ev, s = Base.Checked.sub_with_overflow(a.ev, b.ev)
    isempty(s), FixedMonom(ev, a.deg - b.deg)
end

function monom_is_equal(a::FixedMonom{N}, b::FixedMonom{N}) where {N}
    a.ev == b.ev
end

function monom_create_divmask(
    a::FixedMonom,
    ::Type{Mask},
    ndivvars::Int,
    divmap::Vector{U},
    ndivbits::Int,
    strategy::Symbol
) where {Mask, U}
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
