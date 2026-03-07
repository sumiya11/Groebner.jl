using SmallCollections: FixedVector, SmallVector, fixedvector, sum_fast, dot_fast, bits

struct NibbleNoDeg{N}
    ev::FixedVector{N,UInt8}
end

Base.@propagate_inbounds function monom_exponent(a::NibbleNoDeg, i::Int)
    isodd(i) ? a.ev[i ÷ 2 + 1] & 0x0f : a.ev[i ÷ 2] >> 4
end

lower(a::NibbleNoDeg{N}) where N = a.ev .& 0x0f
upper(a::NibbleNoDeg{N}) where N = a.ev .>> 4
upper_raw(a::NibbleNoDeg{N}) where N = a.ev .& 0xf0

monom_max_vars(a::NibbleNoDeg) = monom_max_vars(typeof(a))

monom_max_vars(::Type{<:NibbleNoDeg{N}}) where N = 2*N
function monom_totaldeg(a::NibbleNoDeg)
    d = reduce(+, lower(a) + upper(a); init = zero(UInt16))
    # d > typemax(UInt8) && __throw_monom_overflow_error()
    d % UInt
end
monom_copy(a::NibbleNoDeg) = a
monom_entrytype(::NibbleNoDeg) = UInt8
monom_entrytype(::Type{<:NibbleNoDeg}) = UInt8  # TODO: shall we indicate somehow that we only have 4 bits per exponent?

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{NibbleNoDeg{N}}, ::Integer) where N
    rand(rng, FixedVector{N,MonomHash})
end

function monom_construct_from_vector(::Type{NibbleNoDeg{N}}, ev::AbstractVector) where N
    all(<=(15), ev) || __throw_monom_overflow_error()
    v = fixedvector(SmallVector{2*N,UInt8}(ev))
    ev = FixedVector{N}(view(v, 1:2:2*N-1)) .| (FixedVector{N}(view(v, 2:2:2*N)) .<< 4)
    NibbleNoDeg(ev)
end

function monom_construct_const(::Type{NibbleNoDeg{N}}, ::Integer) where N
    NibbleNoDeg(zero(FixedVector{N,UInt8}))
end

function monom_hash(a::NibbleNoDeg{N}, b::FixedVector{N}) where N
    dot_fast(a.ev, b)::MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::NibbleNoDeg{N}) where N
    for i in 1:length(tmp)
        @inbounds tmp[i] = monom_exponent(a, i)
    end
    tmp
end

function monom_is_supported_ordering(::Type{<:NibbleNoDeg}, ::Ord) where Ord
    Ord <: Union{DegLex{true}, DegRevLex{true}, InputOrdering}
end

function monom_isless(a::NibbleNoDeg{N}, b::NibbleNoDeg{N}, ::DegLex{true}) where N
    monom_totaldeg(a) != monom_totaldeg(b) && return monom_totaldeg(a) < monom_totaldeg(b)
    for (k, l) in zip(a.ev, b.ev)
        k != l && return bitrotate(k, 4) < bitrotate(l, 4)
    end
    false
end

function monom_isless(a::NibbleNoDeg{N}, b::NibbleNoDeg{N}, ::DegRevLex{true}) where N
    monom_totaldeg(a) < monom_totaldeg(b) || (monom_totaldeg(a) == monom_totaldeg(b) && bits(a.ev) > bits(b.ev))
end

function monom_lcm!(_, a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    NibbleNoDeg(max.(lower(a), lower(b)) .| max.(upper_raw(a), upper_raw(b)))
end

function monom_is_gcd_const(a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    iszero(min.(lower(a), lower(b)) .| min.(upper_raw(a), upper_raw(b)))
end

function monom_product!(_, a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    ev, s = Base.Checked.add_with_overflow(a.ev, b.ev)
    lower_check = (ev .⊻ a.ev .⊻ b.ev) .& 0x10
    if isempty(s) & iszero(lower_check)
        # d = monom_totaldeg(a) + monom_totaldeg(b)
        # d < monom_totaldeg(a) && __throw_monom_overflow_error()
        NibbleNoDeg(ev)
    else
        __throw_monom_overflow_error()
    end
end

function monom_division!(_, a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    NibbleNoDeg(a.ev - b.ev)
end

function monom_is_divisible(a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    first(monom_is_divisible!(a, a, b))  # the first `a` is a dummy
end

function monom_is_divisible!(_, a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    ev, s = Base.Checked.sub_with_overflow(a.ev, b.ev)
    lower_check = iszero((ev .⊻ a.ev .⊻ b.ev) .& 0x10)
    isempty(s) & lower_check, NibbleNoDeg(ev)
end

function monom_is_equal(a::NibbleNoDeg{N}, b::NibbleNoDeg{N}) where N
    a.ev == b.ev
end

function monom_create_divmask(
    a::NibbleNoDeg{N},
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
