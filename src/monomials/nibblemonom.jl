using SmallCollections: FixedVector, SmallVector, fixedvector, sum_fast, dot_fast, bits

struct NibbleMonom{N}
    ev::FixedVector{N,UInt8}
    deg::UInt8
end

function NibbleMonom(ev::FixedVector{N}) where N
    a = NibbleMonom(ev, zero(UInt8))
    d = reduce(+, lower(a) + upper(a); init = zero(UInt16))
    d > typemax(UInt8) && __throw_monom_overflow_error()
    NibbleMonom(ev, d % UInt8)
end

Base.@propagate_inbounds function monom_exponent(a::NibbleMonom, i::Int)
    isodd(i) ? a.ev[i ÷ 2 + 1] & 0x0f : a.ev[i ÷ 2] >> 4
end

lower(a::NibbleMonom{N}) where N = a.ev .& 0x0f
upper(a::NibbleMonom{N}) where N = a.ev .>> 4
upper_raw(a::NibbleMonom{N}) where N = a.ev .& 0xf0

monom_max_vars(a::NibbleMonom) = monom_max_vars(typeof(a))

monom_max_vars(::Type{<:NibbleMonom{N}}) where N = 2*N
monom_totaldeg(a::NibbleMonom) = a.deg % UInt
monom_copy(a::NibbleMonom) = a
monom_entrytype(::NibbleMonom) = UInt8
monom_entrytype(::Type{<:NibbleMonom}) = UInt8  # TODO: shall we indicate somehow that we only have 4 bits per exponent?

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{NibbleMonom{N}}, ::Integer) where N
    rand(rng, FixedVector{N,MonomHash})
end

function monom_construct_from_vector(::Type{NibbleMonom{N}}, ev::AbstractVector) where N
    all(<=(15), ev) || __throw_monom_overflow_error()
    v = fixedvector(SmallVector{2*N,UInt8}(ev))
    ev = FixedVector{N}(view(v, 1:2:2*N-1)) .| (FixedVector{N}(view(v, 2:2:2*N)) .<< 4)
    NibbleMonom(ev)
end

function monom_construct_const(::Type{NibbleMonom{N}}, ::Integer) where N
    NibbleMonom(zero(FixedVector{N,UInt8}))
end

function monom_hash(a::NibbleMonom{N}, b::FixedVector{N}) where N
    dot_fast(a.ev, b)::MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::NibbleMonom{N}) where N
    for i in 1:length(tmp)
        @inbounds tmp[i] = monom_exponent(a, i)
    end
    tmp
end

function monom_is_supported_ordering(::Type{<:NibbleMonom}, ::Ord) where Ord
    Ord <: Union{DegLex{true}, DegRevLex{true}, InputOrdering}
end

function monom_isless(a::NibbleMonom{N}, b::NibbleMonom{N}, ::DegLex{true}) where N
    a.deg != b.deg && return a.deg < b.deg
    for (k, l) in zip(a.ev, b.ev)
        k != l && return bitrotate(k, 4) < bitrotate(l, 4)
    end
    false
end

function monom_isless(a::NibbleMonom{N}, b::NibbleMonom{N}, ::DegRevLex{true}) where N
    a.deg < b.deg || (a.deg == b.deg && bits(a.ev) > bits(b.ev))
end

function monom_lcm!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    NibbleMonom(max.(lower(a), lower(b)) .| max.(upper_raw(a), upper_raw(b)))
end

function monom_is_gcd_const(a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    iszero(min.(lower(a), lower(b)) .| min.(upper_raw(a), upper_raw(b)))
end

function monom_product!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    ev, s = Base.Checked.add_with_overflow(a.ev, b.ev)
    lower_check = (ev .⊻ a.ev .⊻ b.ev) .& 0x10
    if isempty(s) & iszero(lower_check)
        d = a.deg + b.deg
        d < a.deg && __throw_monom_overflow_error()
        NibbleMonom(ev, d)
    else
        __throw_monom_overflow_error()
    end
end

function monom_division!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    NibbleMonom(a.ev - b.ev, a.deg - b.deg)
end

function monom_is_divisible(a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    first(monom_is_divisible!(a, a, b))  # the first `a` is a dummy
end

function monom_is_divisible!(_, a::NibbleMonom{N}, b::NibbleMonom{N}) where N
    ev, s = Base.Checked.sub_with_overflow(a.ev, b.ev)
    lower_check = iszero((ev .⊻ a.ev .⊻ b.ev) .& 0x10)
    isempty(s) & lower_check, NibbleMonom(ev, a.deg - b.deg)
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
