# Defines PackedPairI types for I = 1,2...

# PackedPairI{T, B} implements the interface of a vector
# of small nonegative integers with O(1) vector sum queries. 
# 
# PackedPairI{T, B} is implemented as a sequence of I integer numbers,
# each integer {T} packs a fixed amount of small numbers {B} together.
#
# PackedPairI always stores the sum of vector explicitly.
# PackedPairI always stores in degree-reverse-lex favorable order.

abstract type AbstractPackedPair{T<:Unsigned, B<:Unsigned} end

struct PackedPair1{T<:Unsigned, B<:Unsigned} <: AbstractPackedPair{T, B}
    a1::T    
end

struct PackedPair2{T<:Unsigned, B<:Unsigned} <: AbstractPackedPair{T, B}
    a1::T
    a2::T
end

struct PackedPair3{T<:Unsigned, B<:Unsigned} <: AbstractPackedPair{T, B}
    a1::T
    a2::T
    a3::T
end

# a `p` object can store capacity(p) integers at max. 
capacity(p::AbstractPackedPair) = capacity(typeof(p))

# checks that there is no risk of overflow for `e`.
# If overflow if probable, throws.
function _overflow_check(e::AbstractPackedPair{T, B}) where {T, B}
    _overflow_check(totaldeg(e), B)
end

const _defined_packed_pairs = (
    (:PackedPair1, 1),
    (:PackedPair2, 2),
    (:PackedPair3, 3)
)

# for each PackedPairI define something..
for (op, n) in _defined_packed_pairs
    # define add-on constructors
    @eval begin 
        $op(ev) = $op{UInt64}(ev)
        $op{T}(ev) where {T<:Unsigned} = $op(T, UInt8)(ev)
    end

    # extend `Base.eltype` 
    @eval begin
        Base.eltype(::$op{T, B}) where {T, B} = MonomHash
        Base.eltype(::Type{$op{T, B}}) where {T, B} = MonomHash
    end

    # define `capacity`
    @eval begin
        capacity(::Type{$op{T, B}}) where {T, B} = $n*div(sizeof(T), sizeof(B)) - 1
    end

    # define `totaldeg`
    @eval begin
        @generated function totaldeg(pv::$op{T, B}) where {T, B}
            d = 8*(sizeof(T) - sizeof(B))
            :(pv.a1 >> $d)
        end
    end

    # define `make_hasher`
    @eval begin
        function make_hasher(::Type{$op{T, B}}, n::Integer) where {T, B}
            rand(MonomHash, $n*div(sizeof(T), sizeof(B)))
        end 
    end
end

Base.copy(pv::PackedPair1{T, B}) where {T, B} = PackedPair1{T,B}(pv.a1)
Base.copy(pv::PackedPair2{T, B}) where {T, B} = PackedPair2{T,B}(pv.a1, pv.a2)
Base.copy(pv::PackedPair3{T, B}) where {T, B} = PackedPair3{T,B}(pv.a1, pv.a2, pv.a3)

# Creates an exponent vector of the given type from regular vector `ev`
function make_ev(::Type{PackedPair1{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < epc
    indent = sizeof(T) - degsize(T, B, n)
    a1 = zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _overflow_check(ev[i], B)
        d = T(ev[i])
        a1 = a1 << (sizeof(B)*8)
        a1 = a1 | d
        _overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent*8)
    PackedPair1{T, B}(a1)
end
function make_ev(::Type{PackedPair2{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < 2*epc
    if n < epc
        small = make_ev(PackedPair1{T, B}, ev)
        return PackedPair2{T, B}(small.a1, zero(T))
    end
    indent = sizeof(T) - degsize(T, B, n)
    a1,a2 = zero(T), zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _overflow_check(ev[i], B)
        d = T(ev[i])
        if div(i - 1, epc) == 1
            a1 = a1 << (sizeof(B)*8)
            a1 = a1 | d
        else
            a2 = a2 << (sizeof(B)*8)
            a2 = a2 | d
        end
        _overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent*8)
    PackedPair2{T, B}(a1, a2)
end
function make_ev(::Type{PackedPair3{T, B}}, ev::Vector{U}) where {T, B, U}
    n = length(ev)
    epc = elperchunk(T, B)
    @assert n < 3*epc
    if n < 2*epc
        small = make_ev(PackedPair2{T, B}, ev)
        return PackedPair3{T, B}(small.a1, small.a2, zero(T))
    end
    indent = sizeof(T) - degsize(T, B, n)
    a1,a2,a3 = zero(T), zero(T), zero(T)
    s = zero(T)
    @inbounds for i in n:-1:1
        _overflow_check(ev[i], B)
        d = T(ev[i])
        if div(i - 1, epc) == 2
            a1 = a1 << (sizeof(B)*8)
            a1 = a1 | d
        elseif div(i - 1, epc) == 1
            a2 = a2 << (sizeof(B)*8)
            a2 = a2 | d
        else
            a3 = a3 << (sizeof(B)*8)
            a3 = a3 | d
        end
        _overflow_check(s, B)
        s += d
    end
    a1 |= s << (indent*8)
    PackedPair3{T, B}(a1, a2, a3)
end

# Creates a zero exponent vector of the given type of length n
function make_zero_ev(::Type{PackedPair1{T, B}}, n::Integer) where {T, B}
    PackedPair1{T, B}(zero(T))
end
function make_zero_ev(::Type{PackedPair2{T, B}}, n::Integer) where {T, B}
    PackedPair2{T, B}(zero(T), zero(T))
end
function make_zero_ev(::Type{PackedPair3{T, B}}, n::Integer) where {T, B}
    PackedPair3{T, B}(zero(T), zero(T), zero(T))
end

# Hash of exponent vector `x`
# Must be linear in x:
# hash(x + y) = hash(x) + hash(x)
function Base.hash(x::PackedPair1{T, B}, b::Vector{MH}) where {T, B, MH}
    h = packeddot(x.a1, b, B, 1)
    mod(h, MonomHash)
end
function Base.hash(x::PackedPair2{T, B}, b::Vector{MH}) where {T, B, MH}
    epc = elperchunk(T, B)
    h = packeddot(x.a2, b, B, 0)
    h = h + packeddot(x.a1, view(b, epc+1:length(b)), B, epc - max(epc - 1, length(b) - epc))
    mod(h, MonomHash)
end
function Base.hash(x::PackedPair3{T, B}, b::Vector{MH}) where {T, B, MH}
    epc = elperchunk(T, B)
    h = packeddot(x.a3, b, B, 0)
    h = h + packeddot(x.a2, view(b, epc+1:2*epc), B, 0)
    h = h + packeddot(x.a1, view(b, 2*epc+1:length(b)), B, epc - max(epc - 1, length(b) - 2*epc))
    mod(h, MonomHash)
end

# Creates a regular vector from an exponent vector `pv` 
# and writes the answer to `tmp`
function make_dense!(tmp::Vector{I}, pv::PackedPair1{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    indent = epc - min(epc - 1, length(tmp))
    packedunpack!(tmp, pv.a1, B, indent)
    tmp
end
function make_dense!(tmp::Vector{I}, pv::PackedPair2{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    (length(tmp) < epc) && return make_dense!(tmp, PackedPair1{T, B}(pv.a1))
    indent = 0
    packedunpack!(tmp, pv.a2, B, indent)
    indent = epc - min(epc - 1, length(tmp) - epc)
    packedunpack!(view(tmp, epc+1:length(tmp)), pv.a1, B, indent)
    tmp
end
function make_dense!(tmp::Vector{I}, pv::PackedPair3{T, B}) where {I, T, B}
    epc = elperchunk(T, B)
    (length(tmp) < 2*epc) && return make_dense!(tmp, PackedPair2{T, B}(pv.a1, pv.a2))
    indent = 0
    packedunpack!(tmp, pv.a3, B, indent)
    indent = 0
    packedunpack!(view(tmp, epc+1:2*epc), pv.a2, B, indent)
    indent = epc - min(epc - 1, length(tmp) - 2*epc)
    packedunpack!(view(tmp, 2*epc+1:length(tmp)), pv.a1, B, indent)
    tmp
end

#------------------------------------------------------------------------------
# Monomial orderings implementations 
# for the `PackedPairI` monomial implementation.
# See monoms/orderings.jl for details.

function monom_isless(
        ea::PackedPair1{T, B}, eb::PackedPair1{T, B}, 
        ord::Ord) where {T, B, Ord<:AbstractMonomialOrdering}
    s = div(sizeof(T), sizeof(B))-1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = make_ev(PowerVector{T}, make_dense!(tmp1, ea))
    b = make_ev(PowerVector{T}, make_dense!(tmp2, eb))
    monom_isless(a, b, ord)
end
function monom_isless(
        ea::PackedPair2{T, B}, eb::PackedPair2{T, B}, 
        ord::Ord) where {T, B, Ord<:AbstractMonomialOrdering}    
    s = 2*div(sizeof(T), sizeof(B))-1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = make_ev(PowerVector{T}, make_dense!(tmp1, ea))
    b = make_ev(PowerVector{T}, make_dense!(tmp2, eb))
    monom_isless(a, b, ord)
end
function monom_isless(
        ea::PackedPair3{T, B}, eb::PackedPair3{T, B}, 
        ord::Ord) where {T, B, Ord<:AbstractMonomialOrdering}   
    s = 3*div(sizeof(T), sizeof(B))-1
    tmp1, tmp2 = Vector{T}(undef, s), Vector{T}(undef, s)
    a = make_ev(PowerVector{T}, make_dense!(tmp1, ea))
    b = make_ev(PowerVector{T}, make_dense!(tmp2, eb))
    monom_isless(a, b, ord)
end

function monom_isless(ea::PackedPair1{T, B}, eb::PackedPair1{T, B}, ::DegRevLex) where {T, B}
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

function monom_isless(ea::PackedPair2{T, B}, eb::PackedPair2{T, B}, ::DegRevLex) where {T, B}
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

function monom_isless(ea::PackedPair3{T, B}, eb::PackedPair3{T, B}, ::DegRevLex) where {T, B}
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

#------------------------------------------------------------------------------
# Monomial-Monomial arithmetic.

function monom_lcm!(ec::PackedPair1{T, B}, ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    x, si = packedmax(ea.a1, eb.a1, B, Val(1))
    x += si << ((sizeof(T) - sizeof(B))*8)
    ans = PackedPair1{T, B}(x)
    _overflow_check(ans)
    ans
end
function monom_lcm!(ec::PackedPair2{T, B}, ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    x1, si1 = packedmax(ea.a1, eb.a1, B, Val(1))
    x2, si2 = packedmax(ea.a2, eb.a2, B, Val(0))
    x1 = x1 + ((si1 + si2) << ((sizeof(T) - sizeof(B))*8))
    ans = PackedPair2{T, B}(x1, x2)
    _overflow_check(ans)
    ans
end
function monom_lcm!(ec::PackedPair3{T, B}, ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
    x1, si1 = packedmax(ea.a1, eb.a1, B, Val(1))
    x2, si2 = packedmax(ea.a2, eb.a2, B, Val(0))
    x3, si3 = packedmax(ea.a3, eb.a3, B, Val(0))
    x1 = x1 + ((si1 + si2 + si3) << ((sizeof(T) - sizeof(B))*8))
    ans = PackedPair3{T, B}(x1, x2, x3)
    _overflow_check(ans)
    ans
end

function is_gcd_const(ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    if !packedorth(ea.a1, eb.a1, B, Val(1))
        return false
    end
    return true
end
function is_gcd_const(ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    if !packedorth(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedorth(ea.a2, eb.a2, B, Val(0))
        return false
    end
    return true
end
function is_gcd_const(ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
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

function monom_product!(ec::PackedPair1{T, B}, ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    x = ea.a1 + eb.a1
    ans = PackedPair1{T, B}(x)
    _overflow_check(ans)
    ans
end
function monom_product!(ec::PackedPair2{T, B}, ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    x1 = ea.a1 + eb.a1
    x2 = ea.a2 + eb.a2
    ans = PackedPair2{T, B}(x1, x2)
    _overflow_check(ans)
    ans
end
function monom_product!(ec::PackedPair3{T, B}, ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
    x1 = ea.a1 + eb.a1
    x2 = ea.a2 + eb.a2
    x3 = ea.a3 + eb.a3
    ans = PackedPair3{T, B}(x1, x2, x3)
    _overflow_check(ans)
    ans
end

function monom_division!(ec::PackedPair1{T, B}, ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    x = ea.a1 - eb.a1
    ans = PackedPair1{T, B}(x)
    ans
end
function monom_division!(ec::PackedPair2{T, B}, ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    x1 = ea.a1 - eb.a1
    x2 = ea.a2 - eb.a2
    ans = PackedPair2{T, B}(x1, x2)
    ans
end
function monom_division!(ec::PackedPair3{T, B}, ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
    x1 = ea.a1 - eb.a1
    x2 = ea.a2 - eb.a2
    x3 = ea.a3 - eb.a3
    ans = PackedPair3{T, B}(x1, x2, x3)
    ans
end

function is_monom_divisible(ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    if !packedge(ea.a1, eb.a1, B, Val(1))
        return false
    end
    return true
end
function is_monom_divisible(ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    if !packedge(ea.a1, eb.a1, B, Val(1))
        return false
    end
    if !packedge(ea.a2, eb.a2, B, Val(0))
        return false
    end
    return true
end
function is_monom_divisible(ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
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

function is_monom_divisible!(ec::PackedPair1{T, B}, ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e 
end
function is_monom_divisible!(ec::PackedPair2{T, B}, ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e 
end
function is_monom_divisible!(ec::PackedPair3{T, B}, ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
    ans = is_monom_divisible(ea, eb)
    e = ec
    ans && (e = monom_division!(ec, ea, eb))
    ans, e 
end

function is_monom_elementwise_eq(ea::PackedPair1{T, B}, eb::PackedPair1{T, B}) where {T, B}
    ea.a1 == eb.a1
end
function is_monom_elementwise_eq(ea::PackedPair2{T, B}, eb::PackedPair2{T, B}) where {T, B}
    ea.a1 == eb.a1 && ea.a2 == eb.a2
end
function is_monom_elementwise_eq(ea::PackedPair3{T, B}, eb::PackedPair3{T, B}) where {T, B}
    ea.a1 == eb.a1 && ea.a2 == eb.a2 && ea.a3 == eb.a3
end

#------------------------------------------------------------------------------
# Monomial division masks.
# See f4/hashtable.jl for details.

function monom_divmask(
        e::PackedPair1{T, B},
        DM::Type{Mask},
        ndivvars, divmap,
        ndivbits) where {T, B, Mask}

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    a1 = e.a1
    for i in 1:ndivvars
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B)*8)
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
        e::PackedPair2{T, B},
        DM::Type{Mask},
        ndivvars, divmap,
        ndivbits) where {T, B, Mask}
            
    epc = div(sizeof(T), sizeof(B))

    if ndivvars < epc
        return monom_divmask(PackedPair1{T, B}(e.a1), DM, ndivvars, divmap, ndivbits)
    end

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    
    a2 = e.a2
    for i in 1:epc
        ei = mod(a2, B)
        a2 = a2 >> (sizeof(B)*8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a1 = e.a1
    for i in epc+1:min(2*epc - 1, ndivvars)
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B)*8)
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
    e::PackedPair3{T, B},
    DM::Type{Mask},
    ndivvars, divmap,
    ndivbits) where {T, B, Mask}
    
    epc = elperchunk(T, B)

    if ndivvars < 2*epc
        return monom_divmask(PackedPair2{T, B}(e.a1, e.a2), DM, ndivvars, divmap, ndivbits)
    end

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)

    a3 = e.a3
    for i in 1:epc
        ei = mod(a3, B)
        a3 = a3 >> (sizeof(B)*8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a2 = e.a2
    for i in epc+1:2*epc
        ei = mod(a2, B)
        a2 = a2 >> (sizeof(B)*8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    a1 = e.a1
    for i in 2*epc+1:min(3*epc - 1, ndivvars)
        ei = mod(a1, B)
        a1 = a1 >> (sizeof(B)*8)
        for j in 1:ndivbits
            @inbounds if ei >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end
    
    res
end
