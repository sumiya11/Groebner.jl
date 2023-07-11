# Defines PackedTupleI types for I = 1,2,...,N

# PackedTupleI{T, B} implements the interface of a vector
# of small nonegative integers with O(1) vector sum queries. 
# 
# PackedTupleI{T, B} is implemented as a sequence of I integer numbers,
# each integer {T} packs a fixed amount of small numbers {B} together.
#
# PackedTupleI always stores the sum of vector explicitly.
# PackedTupleI always stores in degree-reverse-lex favorable order.

abstract type AbstractPackedTuple{T <: Unsigned, B <: Unsigned} end

# Checks whether AbstractPackedTuple{T, B} provides a comparator function implementation
# for the given monomial ordering of type `O`
function is_supported_ordering(::Type{APP}, ::O) where {APP <: AbstractPackedTuple, O}
    O <: Union{Lex, DegLex, DegRevLex, InputOrdering}
end

struct PackedTuple1{T <: Unsigned, B <: Unsigned} <: AbstractPackedTuple{T, B}
    a1::T
end

struct PackedTuple2{T <: Unsigned, B <: Unsigned} <: AbstractPackedTuple{T, B}
    a1::T
    a2::T
end

struct PackedTuple3{T <: Unsigned, B <: Unsigned} <: AbstractPackedTuple{T, B}
    a1::T
    a2::T
    a3::T
end

# a `p` object can store max_vars_in_monom(p) integers at max. 
max_vars_in_monom(p::AbstractPackedTuple) = max_vars_in_monom(typeof(p))

# checks that there is no risk of overflow for `e`.
# If overflow if probable, throws.
function _monom_overflow_check(e::AbstractPackedTuple{T, B}) where {T, B}
    _monom_overflow_check(totaldeg(e), B)
end

const _defined_packed_pairs = ((:PackedTuple1, 1), (:PackedTuple2, 2), (:PackedTuple3, 3))

# for each PackedTupleI define something..
for (op, n) in _defined_packed_pairs
    # define add-on constructors
    @eval begin
        $op(ev) = $op{UInt64}(ev)
        $op{T}(ev) where {T <: Unsigned} = $op(T, UInt8)(ev)
    end

    # extend `Base.eltype` 
    @eval begin
        Base.eltype(::$op{T, B}) where {T, B} = MonomHash
        Base.eltype(::Type{$op{T, B}}) where {T, B} = MonomHash
    end

    # define `max_vars_in_monom`
    @eval begin
        max_vars_in_monom(::Type{$op{T, B}}) where {T, B} = $n * div(sizeof(T), sizeof(B)) - 1
    end

    # define `totaldeg`
    @eval begin
        @generated function totaldeg(pv::$op{T, B}) where {T, B}
            d = 8 * (sizeof(T) - sizeof(B))
            :(pv.a1 >> $d)
        end
    end

    # define `construct_hash_vector`
    @eval begin
        function construct_hash_vector(::Type{$op{T, B}}, n::Integer) where {T, B}
            rand(MonomHash, $n * div(sizeof(T), sizeof(B)))
        end
    end
end

Base.copy(pv::PackedTuple1{T, B}) where {T, B} = PackedTuple1{T, B}(pv.a1)
Base.copy(pv::PackedTuple2{T, B}) where {T, B} = PackedTuple2{T, B}(pv.a1, pv.a2)
Base.copy(pv::PackedTuple3{T, B}) where {T, B} = PackedTuple3{T, B}(pv.a1, pv.a2, pv.a3)

# Creates an exponent vector of the given type from regular vector `ev`
function construct_monom(::Type{PackedTuple1{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < epc
    indent = sizeof(T) - degsize(T, B, n)
    a1 = zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _monom_overflow_check(ev[i], B)
        d = T(ev[i])
        a1 = a1 << (sizeof(B) * 8)
        a1 = a1 | d
        _monom_overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent * 8)
    PackedTuple1{T, B}(a1)
end
function construct_monom(::Type{PackedTuple2{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < 2 * epc
    if n < epc
        small = construct_monom(PackedTuple1{T, B}, ev)
        return PackedTuple2{T, B}(small.a1, zero(T))
    end
    indent = sizeof(T) - degsize(T, B, n)
    a1, a2 = zero(T), zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _monom_overflow_check(ev[i], B)
        d = T(ev[i])
        if div(i - 1, epc) == 1
            a1 = a1 << (sizeof(B) * 8)
            a1 = a1 | d
        else
            a2 = a2 << (sizeof(B) * 8)
            a2 = a2 | d
        end
        _monom_overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent * 8)
    PackedTuple2{T, B}(a1, a2)
end
function construct_monom(::Type{PackedTuple3{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < 3 * epc
    if n < 2 * epc
        small = construct_monom(PackedTuple2{T, B}, ev)
        return PackedTuple3{T, B}(small.a1, small.a2, zero(T))
    end
    indent = sizeof(T) - degsize(T, B, n)
    a1, a2, a3 = zero(T), zero(T), zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _monom_overflow_check(ev[i], B)
        d = T(ev[i])
        if div(i - 1, epc) == 2
            a1 = a1 << (sizeof(B) * 8)
            a1 = a1 | d
        elseif div(i - 1, epc) == 1
            a2 = a2 << (sizeof(B) * 8)
            a2 = a2 | d
        else
            a3 = a3 << (sizeof(B) * 8)
            a3 = a3 | d
        end
        _monom_overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent * 8)
    PackedTuple3{T, B}(a1, a2, a3)
end
function make_ev(::Type{PackedPair4{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < 4 * epc
    if n < 3 * epc
        small = make_ev(PackedPair3{T, B}, ev)
        return PackedPair4{T, B}(small.a1, small.a2, small.a3, zero(T))
    end
    indent = sizeof(T) - degsize(T, B, n)
    a1, a2, a3, a4 = zero(T), zero(T), zero(T), zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _overflow_check(ev[i], B)
        d = T(ev[i])
        if div(i - 1, epc) == 3
            a1 = a1 << (sizeof(B) * 8)
            a1 = a1 | d
        elseif div(i - 1, epc) == 2
            a2 = a2 << (sizeof(B) * 8)
            a2 = a2 | d
        elseif div(i - 1, epc) == 1
            a3 = a3 << (sizeof(B) * 8)
            a3 = a3 | d
        else
            a4 = a4 << (sizeof(B) * 8)
            a4 = a4 | d
        end
        _overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent * 8)
    PackedPair4{T, B}(a1, a2, a3, a4)
end

# Creates a zero exponent vector of the given type of length n
function construct_const_monom(::Type{PackedTuple1{T, B}}, n::Integer) where {T, B}
    PackedTuple1{T, B}(zero(T))
end
function construct_const_monom(::Type{PackedTuple2{T, B}}, n::Integer) where {T, B}
    PackedTuple2{T, B}(zero(T), zero(T))
end
function construct_const_monom(::Type{PackedTuple3{T, B}}, n::Integer) where {T, B}
    PackedTuple3{T, B}(zero(T), zero(T), zero(T))
end
function make_zero_ev(::Type{PackedPair4{T, B}}, n::Integer) where {T, B}
    PackedPair4{T, B}(zero(T), zero(T), zero(T), zero(T))
end

# Hash of exponent vector `x`
# Must be linear in x:
# monom_hash(x + y) = monom_hash(x) + monom_hash(x)
function monom_hash(x::PackedTuple1{T, B}, b::Vector{MH}) where {T, B, MH}
    h = packeddot(x.a1, b, B, 1)
    mod(h, MonomHash)
end
function monom_hash(x::PackedTuple2{T, B}, b::Vector{MH}) where {T, B, MH}
    epc = elperchunk(T, B)
    h = packeddot(x.a2, b, B, 0)
    h =
        h + packeddot(
            x.a1,
            view(b, (epc + 1):length(b)),
            B,
            epc - max(epc - 1, length(b) - epc)
        )
    mod(h, MonomHash)
end
function monom_hash(x::PackedTuple3{T, B}, b::Vector{MH}) where {T, B, MH}
    epc = elperchunk(T, B)
    h = packeddot(x.a3, b, B, 0)
    h = h + packeddot(x.a2, view(b, (epc + 1):(2 * epc)), B, 0)
    h =
        h + packeddot(
            x.a1,
            view(b, (2 * epc + 1):length(b)),
            B,
            epc - max(epc - 1, length(b) - 2 * epc)
        )
    mod(h, MonomHash)
end
function Base.hash(x::PackedPair4{T, B}, b::Vector{MH}) where {T, B, MH}
    epc = elperchunk(T, B)
    h = packeddot(x.a4, b, B, 0)
    h = h + packeddot(x.a3, view(b, (epc + 1):(2 * epc)), B, 0)
    h = h + packeddot(x.a2, view(b, (2epc + 1):(3 * epc)), B, 0)
    h =
        h + packeddot(
            x.a1,
            view(b, (3 * epc + 1):length(b)),
            B,
            epc - max(epc - 1, length(b) - 3 * epc)
        )
    mod(h, MonomHash)
end

# Creates a regular vector from an exponent vector `pv` 
# and writes the answer to `tmp`
function monom_to_dense_vector!(tmp::Vector{I}, pv::PackedTuple1{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    indent = epc - min(epc - 1, length(tmp))
    packedunpack!(tmp, pv.a1, B, indent)
    tmp
end
function monom_to_dense_vector!(tmp::Vector{I}, pv::PackedTuple2{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    (length(tmp) < epc) && return monom_to_dense_vector!(tmp, PackedTuple1{T, B}(pv.a1))
    indent = 0
    packedunpack!(tmp, pv.a2, B, indent)
    indent = epc - min(epc - 1, length(tmp) - epc)
    packedunpack!(view(tmp, (epc + 1):length(tmp)), pv.a1, B, indent)
    tmp
end
function monom_to_dense_vector!(tmp::Vector{I}, pv::PackedTuple3{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    (length(tmp) < 2 * epc) && return monom_to_dense_vector!(tmp, PackedTuple2{T, B}(pv.a1, pv.a2))
    indent = 0
    packedunpack!(tmp, pv.a3, B, indent)
    indent = 0
    packedunpack!(view(tmp, (epc + 1):(2 * epc)), pv.a2, B, indent)
    indent = epc - min(epc - 1, length(tmp) - 2 * epc)
    packedunpack!(view(tmp, (2 * epc + 1):length(tmp)), pv.a1, B, indent)
    tmp
end
function make_dense!(tmp::Vector{I}, pv::PackedPair4{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    (length(tmp) < 3 * epc) &&
        return make_dense!(tmp, PackedPair3{T, B}(pv.a1, pv.a2, pv.a3))
    indent = 0
    packedunpack!(tmp, pv.a4, B, indent)
    indent = 0
    packedunpack!(view(tmp, (epc + 1):(2 * epc)), pv.a3, B, indent)
    indent = 0
    packedunpack!(view(tmp, (2 * epc + 1):(3 * epc)), pv.a2, B, indent)
    indent = epc - min(epc - 1, length(tmp) - 3 * epc)
    packedunpack!(view(tmp, (3 * epc + 1):length(tmp)), pv.a1, B, indent)
    tmp
end

#------------------------------------------------------------------------------
# Monomial orderings implementations 
# for the `PackedTupleI` monomial implementation.
# See monoms/orderings.jl for details.

function monom_isless(
    ea::PackedTuple1{T, B},
    eb::PackedTuple1{T, B},
    ord::Ord
) where {T, B, Ord <: AbstractMonomialOrdering}
    s = div(sizeof(T), sizeof(B)) - 1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = construct_monom(ExponentVector{T}, monom_to_dense_vector!(tmp1, ea))
    b = construct_monom(ExponentVector{T}, monom_to_dense_vector!(tmp2, eb))
    monom_isless(a, b, ord)
end
function monom_isless(
    ea::PackedTuple2{T, B},
    eb::PackedTuple2{T, B},
    ord::Ord
) where {T, B, Ord <: AbstractMonomialOrdering}
    s = 2 * div(sizeof(T), sizeof(B)) - 1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = construct_monom(ExponentVector{T}, monom_to_dense_vector!(tmp1, ea))
    b = construct_monom(ExponentVector{T}, monom_to_dense_vector!(tmp2, eb))
    monom_isless(a, b, ord)
end
function monom_isless(
    ea::PackedTuple3{T, B},
    eb::PackedTuple3{T, B},
    ord::Ord
) where {T, B, Ord <: AbstractMonomialOrdering}
    s = 3 * div(sizeof(T), sizeof(B)) - 1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = construct_monom(ExponentVector{T}, monom_to_dense_vector!(tmp1, ea))
    b = construct_monom(ExponentVector{T}, monom_to_dense_vector!(tmp2, eb))
    monom_isless(a, b, ord)
end
function monom_isless(
    ea::PackedPair4{T, B},
    eb::PackedPair4{T, B},
    ord::Ord
) where {T, B, Ord <: AbstractMonomialOrdering}
    s = 4 * div(sizeof(T), sizeof(B)) - 1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = make_ev(PowerVector{T}, make_dense!(tmp1, ea))
    b = make_ev(PowerVector{T}, make_dense!(tmp2, eb))
    monom_isless(a, b, ord)
end

function monom_isless(
    ea::PackedTuple1{T, B},
    eb::PackedTuple1{T, B},
    ::DegRevLex
) where {T, B}
    da, db = totaldeg(ea), totaldeg(eb)
    if da < db
        return true
    end
    if da > db
        return false
    end

    if ea.a1 <= eb.a1
        return false
    else
        return true
    end
end

function monom_isless(
    ea::PackedTuple2{T, B},
    eb::PackedTuple2{T, B},
    ::DegRevLex
) where {T, B}
    da, db = totaldeg(ea), totaldeg(eb)
    if da < db
        return true
    end
    if da > db
        return false
    end

    if ea.a1 == eb.a1
        return !(ea.a2 <= eb.a2)
    else
        return !(ea.a1 <= eb.a1)
    end
end

function monom_isless(
    ea::PackedTuple3{T, B},
    eb::PackedTuple3{T, B},
    ::DegRevLex
) where {T, B}
    da, db = totaldeg(ea), totaldeg(eb)
    if da < db
        return true
    end
    if da > db
        return false
    end

    if ea.a1 == eb.a1
        if ea.a2 == eb.a2
            return !(ea.a3 <= eb.a3)
        else
            return !(ea.a2 <= eb.a2)
        end
    else
        return !(ea.a1 <= eb.a1)
    end
end

function monom_isless(
    ea::PackedPair4{T, B},
    eb::PackedPair4{T, B},
    ::DegRevLex
) where {T, B}
    da, db = totaldeg(ea), totaldeg(eb)
    if da < db
        return true
    end
    if da > db
        return false
    end

    if ea.a1 == eb.a1
        if ea.a2 == eb.a2
            if ea.a3 == eb.a3
                return !(ea.a4 <= eb.a4)
            else
                return !(ea.a3 <= eb.a3)
            end
        else
            return !(ea.a2 <= eb.a2)
        end
    else
        return !(ea.a1 <= eb.a1)
    end
end

#------------------------------------------------------------------------------
# Monomial-Monomial arithmetic.

function monom_lcm!(
    ec::PackedTuple1{T, B},
    ea::PackedTuple1{T, B},
    eb::PackedTuple1{T, B}
) where {T, B}
    x, si = packedmax(ea.a1, eb.a1, B, Val(1))
    x += si << ((sizeof(T) - sizeof(B)) * 8)
    ans = PackedTuple1{T, B}(x)
    _monom_overflow_check(ans)
    ans
end
function monom_lcm!(
    ec::PackedTuple2{T, B},
    ea::PackedTuple2{T, B},
    eb::PackedTuple2{T, B}
) where {T, B}
    x1, si1 = packedmax(ea.a1, eb.a1, B, Val(1))
    x2, si2 = packedmax(ea.a2, eb.a2, B, Val(0))
    x1 = x1 + ((si1 + si2) << ((sizeof(T) - sizeof(B)) * 8))
    ans = PackedTuple2{T, B}(x1, x2)
    _monom_overflow_check(ans)
    ans
end
function monom_lcm!(
    ec::PackedTuple3{T, B},
    ea::PackedTuple3{T, B},
    eb::PackedTuple3{T, B}
) where {T, B}
    x1, si1 = packedmax(ea.a1, eb.a1, B, Val(1))
    x2, si2 = packedmax(ea.a2, eb.a2, B, Val(0))
    x3, si3 = packedmax(ea.a3, eb.a3, B, Val(0))
    x1 = x1 + ((si1 + si2 + si3) << ((sizeof(T) - sizeof(B)) * 8))
    ans = PackedTuple3{T, B}(x1, x2, x3)
    _monom_overflow_check(ans)
    ans
end
function monom_lcm!(
    ec::PackedPair4{T, B},
    ea::PackedPair4{T, B},
    eb::PackedPair4{T, B}
) where {T, B}
    x1, si1 = packedmax(ea.a1, eb.a1, B, Val(1))
    x2, si2 = packedmax(ea.a2, eb.a2, B, Val(0))
    x3, si3 = packedmax(ea.a3, eb.a3, B, Val(0))
    x4, si4 = packedmax(ea.a4, eb.a4, B, Val(0))
    x1 = x1 + ((si1 + si2 + si3 + si4) << ((sizeof(T) - sizeof(B)) * 8))
    ans = PackedPair4{T, B}(x1, x2, x3, x4)
    _overflow_check(ans)
    ans
end

function is_gcd_const(ea::PackedTuple1{T, B}, eb::PackedTuple1{T, B}) where {T, B}
    if !packedorth(ea.a1, eb.a1, B, Val(1))
        return false
    end
    return true
end
function is_gcd_const(ea::PackedTuple2{T, B}, eb::PackedTuple2{T, B}) where {T, B}
    if !packedorth(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedorth(ea.a2, eb.a2, B, Val(0))
        return false
    end
    return true
end
function is_gcd_const(ea::PackedTuple3{T, B}, eb::PackedTuple3{T, B}) where {T, B}
    if !packedorth(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedorth(ea.a2, eb.a2, B, Val(0))
        return false
    end
    if !packedorth(ea.a3, eb.a3, B, Val(0))
        return false
    end
    return true
end
function is_gcd_const(ea::PackedPair4{T, B}, eb::PackedPair4{T, B}) where {T, B}
    if !packedorth(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedorth(ea.a2, eb.a2, B, Val(0))
        return false
    end
    if !packedorth(ea.a3, eb.a3, B, Val(0))
        return false
    end
    if !packedorth(ea.a4, eb.a4, B, Val(0))
        return false
    end
    return true
end

function monom_product!(
    ec::PackedTuple1{T, B},
    ea::PackedTuple1{T, B},
    eb::PackedTuple1{T, B}
) where {T, B}
    x = ea.a1 + eb.a1
    ans = PackedTuple1{T, B}(x)
    _monom_overflow_check(ans)
    ans
end
function monom_product!(
    ec::PackedTuple2{T, B},
    ea::PackedTuple2{T, B},
    eb::PackedTuple2{T, B}
) where {T, B}
    x1 = ea.a1 + eb.a1
    x2 = ea.a2 + eb.a2
    ans = PackedTuple2{T, B}(x1, x2)
    _monom_overflow_check(ans)
    ans
end
function monom_product!(
    ec::PackedTuple3{T, B},
    ea::PackedTuple3{T, B},
    eb::PackedTuple3{T, B}
) where {T, B}
    x1 = ea.a1 + eb.a1
    x2 = ea.a2 + eb.a2
    x3 = ea.a3 + eb.a3
    ans = PackedTuple3{T, B}(x1, x2, x3)
    _monom_overflow_check(ans)
    ans
end
function monom_product!(
    ec::PackedPair4{T, B},
    ea::PackedPair4{T, B},
    eb::PackedPair4{T, B}
) where {T, B}
    x1 = ea.a1 + eb.a1
    x2 = ea.a2 + eb.a2
    x3 = ea.a3 + eb.a3
    x4 = ea.a4 + eb.a4
    ans = PackedPair4{T, B}(x1, x2, x3, x4)
    _overflow_check(ans)
    ans
end

function monom_division!(
    ec::PackedTuple1{T, B},
    ea::PackedTuple1{T, B},
    eb::PackedTuple1{T, B}
) where {T, B}
    x = ea.a1 - eb.a1
    ans = PackedTuple1{T, B}(x)
    ans
end
function monom_division!(
    ec::PackedTuple2{T, B},
    ea::PackedTuple2{T, B},
    eb::PackedTuple2{T, B}
) where {T, B}
    x1 = ea.a1 - eb.a1
    x2 = ea.a2 - eb.a2
    ans = PackedTuple2{T, B}(x1, x2)
    ans
end
function monom_division!(
    ec::PackedTuple3{T, B},
    ea::PackedTuple3{T, B},
    eb::PackedTuple3{T, B}
) where {T, B}
    x1 = ea.a1 - eb.a1
    x2 = ea.a2 - eb.a2
    x3 = ea.a3 - eb.a3
    ans = PackedTuple3{T, B}(x1, x2, x3)
    ans
end
function monom_division!(
    ec::PackedPair4{T, B},
    ea::PackedPair4{T, B},
    eb::PackedPair4{T, B}
) where {T, B}
    x1 = ea.a1 - eb.a1
    x2 = ea.a2 - eb.a2
    x3 = ea.a3 - eb.a3
    x4 = ea.a4 - eb.a4
    ans = PackedPair4{T, B}(x1, x2, x3, x4)
    ans
end

function is_monom_divisible(ea::PackedTuple1{T, B}, eb::PackedTuple1{T, B}) where {T, B}
    if !packedge(ea.a1, eb.a1, B, Val(1))
        return false
    end
    return true
end
function is_monom_divisible(ea::PackedTuple2{T, B}, eb::PackedTuple2{T, B}) where {T, B}
    if !packedge(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedge(ea.a2, eb.a2, B, Val(0))
        return false
    end
    return true
end
function is_monom_divisible(ea::PackedTuple3{T, B}, eb::PackedTuple3{T, B}) where {T, B}
    if !packedge(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedge(ea.a2, eb.a2, B, Val(0))
        return false
    end
    if !packedge(ea.a3, eb.a3, B, Val(0))
        return false
    end
    return true
end
function is_monom_divisible(ea::PackedPair4{T, B}, eb::PackedPair4{T, B}) where {T, B}
    if !packedge(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedge(ea.a2, eb.a2, B, Val(0))
        return false
    end
    if !packedge(ea.a3, eb.a3, B, Val(0))
        return false
    end
    if !packedge(ea.a4, eb.a4, B, Val(0))
        return false
    end
    return true
end

function is_monom_divisible!(
    ec::PackedTuple1{T, B},
    ea::PackedTuple1{T, B},
    eb::PackedTuple1{T, B}
) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e
end
function is_monom_divisible!(
    ec::PackedTuple2{T, B},
    ea::PackedTuple2{T, B},
    eb::PackedTuple2{T, B}
) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e
end
function is_monom_divisible!(
    ec::PackedTuple3{T, B},
    ea::PackedTuple3{T, B},
    eb::PackedTuple3{T, B}
) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e
end
function is_monom_divisible!(
    ec::PackedPair4{T, B},
    ea::PackedPair4{T, B},
    eb::PackedPair4{T, B}
) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e
end

function is_monom_elementwise_eq(ea::PackedTuple1{T, B}, eb::PackedTuple1{T, B}) where {T, B}
    ea.a1 == eb.a1
end
function is_monom_elementwise_eq(ea::PackedTuple2{T, B}, eb::PackedTuple2{T, B}) where {T, B}
    ea.a1 == eb.a1 && ea.a2 == eb.a2
end
function is_monom_elementwise_eq(ea::PackedTuple3{T, B}, eb::PackedTuple3{T, B}) where {T, B}
    ea.a1 == eb.a1 && ea.a2 == eb.a2 && ea.a3 == eb.a3
end
function is_monom_elementwise_eq(ea::PackedPair4{T, B}, eb::PackedPair4{T, B}) where {T, B}
    ea.a1 == eb.a1 && ea.a2 == eb.a2 && ea.a3 == eb.a3 && ea.a4 == eb.a4
end

#------------------------------------------------------------------------------
# Monomial division masks.
# See f4/hashtable.jl for details.

function monom_divmask(
    e::PackedTuple1{T, B},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, B, Mask}
    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    a1 = e.a1
    for i in 1:ndivvars
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end

function monom_divmask(
    e::PackedTuple2{T, B},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, B, Mask}
    epc = div(sizeof(T), sizeof(B))

    if ndivvars < epc
        return monom_divmask(PackedTuple1{T, B}(e.a1), DM, ndivvars, divmap, ndivbits)
    end

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)

    a2 = e.a2
    for i in 1:epc
        ei = mod(a2, B)
        a2 = a2 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a1 = e.a1
    for i in (epc + 1):min(2 * epc - 1, ndivvars)
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end

function monom_divmask(
    e::PackedTuple3{T, B},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, B, Mask}
    epc = elperchunk(T, B)

    if ndivvars < 2 * epc
        return monom_divmask(PackedTuple2{T, B}(e.a1, e.a2), DM, ndivvars, divmap, ndivbits)
    end

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)

    a3 = e.a3
    for i in 1:epc
        ei = mod(a3, B)
        a3 = a3 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a2 = e.a2
    for i in (epc + 1):(2 * epc)
        ei = mod(a2, B)
        a2 = a2 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a1 = e.a1
    for i in (2 * epc + 1):min(3 * epc - 1, ndivvars)
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end

function monom_divmask(
    e::PackedPair4{T, B},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, B, Mask}
    epc = elperchunk(T, B)

    if ndivvars < 3 * epc
        return monom_divmask(
            PackedPair3{T, B}(e.a1, e.a2, e.a3),
            DM,
            ndivvars,
            divmap,
            ndivbits
        )
    end

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)

    a4 = e.a4
    for i in 1:epc
        ei = mod(a3, B)
        a4 = a4 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a3 = e.a3
    for i in (epc + 1):(2 * epc)
        ei = mod(a3, B)
        a3 = a3 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a2 = e.a2
    for i in (2 * epc + 1):(3 * epc)
        ei = mod(a2, B)
        a2 = a2 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a1 = e.a1
    for i in (3 * epc + 1):min(4 * epc - 1, ndivvars)
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B) * 8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end
