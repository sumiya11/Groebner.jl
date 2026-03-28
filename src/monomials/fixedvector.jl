using SmallCollections: FixedVector, SmallVector, fixedvector, setindex, sum_fast, bits

monom_max_vars(a::FixedVector) = monom_max_vars(typeof(a))

function monom_overflow_check(a::FixedVector)
    c = ~a[end]
    signbit(signed(c)) && __throw_monom_overflow_error()
end

monom_max_vars(::Type{A}) where A <: FixedVector = capacity(A) - 1
monom_totaldeg(a::FixedVector{N,T}) where {N,T} = ~a[end] % T
monom_copy(a::FixedVector) = a

function monom_construct_hash_vector(rng::AbstractRNG, ::Type{<:FixedVector{N}}, ::Integer) where N
    setindex(rand(rng, FixedVector{N,MonomHash}), zero(MonomHash), N)
end

function monom_construct_from_vector(::Type{FixedVector{N,T}}, ev::AbstractVector) where {N,T}
    d = sum(ev)
    d <= typemax(signed(T)) || __throw_monom_overflow_error()
    a = fixedvector(SmallVector{N,T}(ev))
    setindex(a, ~(d % T), N)
end

function monom_construct_const(::Type{FixedVector{N,T}}, ::Integer) where {N,T}
    setindex(zero(FixedVector{N,T}), ~zero(T), N)
end

function monom_hash(a::FixedVector{N}, b::FixedVector{N}) where N
    dot_fast(a, b) % MonomHash
end

function monom_to_vector!(tmp::AbstractVector, a::FixedVector)
    copyto!(tmp, 1, a, 1, length(tmp))
end

function monom_is_supported_ordering(::Type{<:FixedVector}, ::Ord) where Ord
    Ord <: Union{DegLex{true}, DegRevLex{true}, InputOrdering}
end

function monom_isless(a::FixedVector{N}, b::FixedVector{N}, ::DegLex{true}) where N
    a[end] > b[end] || (a[end] == b[end] && a < b)
end

function monom_isless(a::FixedVector{N}, b::FixedVector{N}, ::DegRevLex{true}) where N
    if sizeof(a) <= 64
        bits(a) > bits(b)
    else
        aev = reinterpret(UInt64, a)
        bev = reinterpret(UInt64, b)
        @inbounds for i in length(aev):-1:1
            aev[i] != bev[i] && return aev[i] > bev[i]
        end
        false
    end
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
    m = bits(a .< b)
    l = one(m) << (N-1)
    iszero(m & ~l)
end

function monom_is_divisible!(_, a::FixedVector{N}, b::FixedVector{N}) where N
    monom_is_divisible(a, b) ? (true, monom_division!(a, a, b)) : (false, a)
    # the first `a` argument for `monom_division!` is a dummy
end

function monom_is_equal(a::FixedVector{N}, b::FixedVector{N}) where N
    a == b
end

function monom_create_divmask(
    e::FixedVector,
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
            if e[i] >= divmap[ctr]
                res |= o << ((ctr - 1) % UInt)
            end
            ctr += o
        end
    end
    res
end
