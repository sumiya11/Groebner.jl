using SmallCollections: FixedVector, SmallVector, fixedvector, setindex, sum_fast, bits, bitsize

monom_max_vars(a::FixedVector) = monom_max_vars(typeof(a))

function monom_overflow_check(a::FixedVector{N,T}) where {N,T}
    c = ~a[end]
    signbit(signed(c)) && __throw_monom_overflow_error(c, T)
end

monom_max_vars(::Type{A}) where A <: FixedVector = capacity(A) - 1
monom_totaldeg(a::FixedVector) = ~a[end] % UInt
monom_copy(a::FixedVector) = a
monom_entrytype(a::FixedVector) = eltype(a)
monom_entrytype(::Type{A}) where A <: FixedVector = eltype(A)

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{FixedVector{N,T}}, ::Integer) where {N,T}
    setindex(rand(FixedVector{N,MonomHash}), zero(MonomHash), N)
end

function monom_construct_from_vector(::Type{FixedVector{N,T}}, ev::AbstractVector) where {N,T}
    a = fixedvector(@inbounds SmallVector{N,T}(ev))
    setindex(a, ~sum_fast(a), N)
end

function monom_construct_const(::Type{FixedVector{N,T}}, ::Integer) where {N,T}
    setindex(zero(FixedVector{N,T}), ~zero(T), N)
end

function monom_hash(a::FixedVector{N}, b::FixedVector{N}) where N
    sum_fast(map(*, a, b)) % MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::FixedVector)
    copyto!(tmp, 1, a, 1, length(tmp))
end

function monom_is_supported_ordering(::Type{<:FixedVector}, ::Ord) where Ord
    Ord <: Union{DegRevLex{true}, InputOrdering}
end

function monom_isless(a::FixedVector{N}, b::FixedVector{N}, ::DegRevLex{true}) where N
    bits(a) > bits(b)
end

function monom_lcm!(_, a::FixedVector{N,T}, b::FixedVector{N,T}) where {N,T}
    u = FixedVector{N,T}(ntuple(i -> i == N ? zero(T) : ~zero(T), Val(N)))
    c = max.(a, b) .& u
    c = setindex(c, ~sum_fast(c), N)
    monom_overflow_check(c)
    c
end

function monom_is_gcd_const(a::FixedVector{N}, b::FixedVector{N}) where N
    a1 = bits(map(!iszero, a))
    b1 = bits(map(!iszero, b))
    l = one(a1) << (N-1)
    iszero(a1 & b1 & ~l)
end

function monom_product!(_, a::FixedVector{N,T}, b::FixedVector{N,T}) where {N,T}
    u = FixedVector{N,T}(ntuple(k -> T(k == N), Val(N)))
    c = a + b + u
    monom_overflow_check(c)
    c
end

function monom_division!(_, a::FixedVector{N,T}, b::FixedVector{N,T}) where {N,T}
    u = FixedVector{N,T}(ntuple(k -> T(k == N), Val(N)))
    c = a - b - u
end

function monom_is_divisible(a::FixedVector{N}, b::FixedVector{N}) where N
    m = bits(map(<, a, b))
    l = one(m) << (N-1)
    iszero(m & ~l)
end

function monom_is_divisible!(c::FixedVector{N}, a::FixedVector{N}, b::FixedVector{N}) where N
    monom_is_divisible(a, b) ? (true, monom_division!(c, a, b)) : (false, c)
end

function monom_is_equal(a::FixedVector{N}, b::FixedVector{N}) where N
    a == b
end

function monom_create_divmask(
    e::FixedVector{N,T},
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
            if e[i + 1] >= divmap[ctr]
                res |= o << ((ctr - 1) % UInt)
            end
            ctr += o
        end
    end
    res
end
