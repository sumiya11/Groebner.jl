using SmallCollections: FixedVector, SmallVector, fixedvector, sum_fast, bits

struct NibbleMonom{N}
    ev::FixedVector{N,UInt8}
    deg::UInt8
end

function NibbleMonom(ev::FixedVector{N}) where N
    a = NibbleMonom(ev, zero(UInt8))
    d = sum_fast(lower(a)) + sum_fast(upper(a))
    a = NibbleMonom(ev, d)
    monom_overflow_check(a)
    a
end

Base.@propagate_inbounds function monom_exponent(a::NibbleMonom, i::Int)
    isodd(i) ? a.ev[i ÷ 2 + 1] & 0x0f : a.ev[i ÷ 2] >> 4
end

lower(a::NibbleMonom{N}) where N = a.ev .& 0x0f
upper(a::NibbleMonom{N}) where N = a.ev .>> 4
upper_raw(a::NibbleMonom{N}) where N = a.ev .& 0xf0

monom_max_vars(a::NibbleMonom) = monom_max_vars(typeof(a))

function monom_overflow_check(a::NibbleMonom)
    signbit(signed(a.deg)) && __throw_monom_overflow_error(a.deg, UInt8)
end

monom_max_vars(::Type{<:NibbleMonom{N}}) where N = 2*N
monom_totaldeg(a::NibbleMonom) = a.deg % UInt
monom_copy(a::NibbleMonom) = a
monom_entrytype(::NibbleMonom) = UInt8
monom_entrytype(::Type{<:NibbleMonom}) = UInt8  # TODO: shall we indicate somehow that we only have 4 bits per exponent?

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{NibbleMonom{N}}, ::Integer) where N
    rand(rng, FixedVector{N,MonomHash})
end

function monom_construct_from_vector(::Type{NibbleMonom{N}}, ev::AbstractVector) where N
    v = fixedvector(SmallVector{2*N,UInt8}(ev))
    ev = FixedVector{N}(view(v, 1:2:2*N-1)) .| (FixedVector{N}(view(v, 2:2:2*N)) .<< 4)
    NibbleMonom(ev)
end

function monom_construct_const(::Type{NibbleMonom{N}}, ::Integer) where N
    NibbleMonom(zero(FixedVector{N,UInt8}))
end

function monom_hash(a::NibbleMonom{N}, b::FixedVector{N}) where N
    sum_fast(map(*, a.ev, b))  # TODO: can we combine two exponents like this?
end

function monom_to_vector!(tmp::AbstractVector, a::NibbleMonom{N}) where N
    #=
    u = FixedVector{N,UInt8}(ntuple(Returns(0x0f), Val(N)))
    copyto!(@view(tmp[1:2:end]), 1, lower(a), 1, (length(tmp)+1) ÷ 2)
    copyto!(@view(tmp[2:2:end]), 1, upper(a), 1, length(tmp) ÷ 2)
    =#
    for i in 1:length(tmp)
        @inbounds tmp[i] = monom_exponent(a, i)
    end
    tmp
end

function monom_is_supported_ordering(::Type{<:NibbleMonom}, ::Ord) where Ord
    Ord <: Union{DegRevLex{true}, InputOrdering}
end

function monom_isless(a::NibbleMonom{N}, b::NibbleMonom{N}, ::DegRevLex{true}) where N
    a.deg < b.deg || (a.deg == b.deg && bits(a.ev) > bits(b.ev))
end

function monom_lcm!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    NibbleMonom(max.(lower(a), lower(b)) .| max.(upper_raw(a), upper_raw(b)))
end

function monom_is_gcd_const(a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    p = bits(map(!iszero, lower(a))) & bits(map(!iszero, lower(b)))
    q = bits(map(!iszero, upper_raw(a))) & bits(map(!iszero, upper_raw(b)))
    iszero(p | q)
end

function monom_product!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    ev = a.ev + b.ev
    iszero((ev .⊻ a.ev .⊻ b.ev) .& 0x10) || error("overflow")  # TODO: check overflow of the high nibble
    c = NibbleMonom(ev, a.deg + b.deg)
    monom_overflow_check(c)
    c
end

function monom_division!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    NibbleMonom(a.ev - b.ev, a.deg - b.deg)
end

function monom_is_divisible(a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    all(map(>=, lower(a), lower(b)) .& map(>=, upper_raw(a), upper_raw(b)))  # TODO: use checked arithmetic?
end

function monom_is_divisible!(c::NibbleMonom{N}, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    monom_is_divisible(a, b) ? (true, monom_division!(c, a, b)) : (false, c)  # TODO: use checked arithmetic
end

function monom_is_equal(a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    a.ev == b.ev
end

function monom_create_divmask(
    a::NibbleMonom{N},
    ::Type{Mask},
    ndivvars::Int,
    divmap::Vector{U},
    ndivbits::Int,
    strategy::Symbol
) where {N, Mask, U}
    @invariant strategy == :first_variables
    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    @inbounds for i in 1:ndivvars
        for _ in 1:ndivbits
            if monom_exponent(a, i) >= divmap[ctr]
                res |= o << ((ctr - 1) % UInt)
            end
            ctr += o
        end
    end
    res
end
