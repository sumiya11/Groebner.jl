# This file implements a new and simple monomial type 
# compatible with Groebner.jl

# This scary line used to load and directly access some 
# internal types in the *current local branch* of Groebner.jl
include((@__DIR__)*"/../../src/Groebner.jl")

#=
 What should be defined here in order for this new monomial type to work ?
=#

#=
1. The monomial type itself should be defined.

   Note that this particular definition suggests
   that objects of type MyNewMonom will be immutable, 
   and the vector `degrees` inside of it is mutable.
=#
struct MyNewMonom
    # a structure with a single vector of monomial degrees,
    # given as a vector of nonnegative integers of type UInt64
    degrees::Vector{UInt64}
end

#=
2. Utility functions should be defined (or, extended).
    These include:
    - _overflow_check(mnm::MyNewMonom)
        this could always return false,
        meaning that there is no risk of exponent overflow in MyNewMonom
    - capacity(::Type{MyNewMonom}),
        the maximum number of variables that can be safely stored in MyNewMonom.
    - capacity(mnm::MyNewMonom),
        -//-
    - totaldeg(mnm::MyNewMonom),
        the total degree of `mnm`
    - powertype(::Type{MyNewMonom})
        this should be either the type of the entry of the monomial (UInt64, in our case),
        or something else (and then, monomial hashing will depend on this type)
    - powertype(p::MyNewMonom)
        -//-
 And
    - make_zero_ev(::Type{MyNewMonom}, n::Integer),
        Creates a zero MyNewMonom object capable of storing degrees of n variables
    - make_ev(::Type{MyNewMonom}, ev::Vector),
        Creates MyNewMonom object with the data from the given vector `ev`
    - make_dense!(tmp::Vector, mnm::MyNewMonom),
        Fills the vector `tmp` with the degrees from the MyNewMonom object `mnm`
    - make_hasher(::Type{MyNewMonom}, n::Integer),
        Creates a hash vector (which should also be of type MyNewMonom),
        which would be used for hashing other MyNewMonom objects
    - Base.hash(x::MyNewMonom, b::MyNewMonom),
        Computes the hash of the MyNewMonom object `x` with the hash vector `b`.
        This function must be linear in `x`.
=#

totaldeg(mnm::MyNewMonom) = sum(mnm.degrees)
_overflow_check(mnm::MyNewMonom) = false

capacity(::Type{MyNewMonom}) = 2^32
capacity(mnm::MyNewMonom) = capacity(typeof(mnm))

powertype(::Type{MyNewMonom}) = UInt64
powertype(::MyNewMonom) = UInt64

make_zero_ev(::Type{MyNewMonom}, n::Integer) = MyNewMonom(zeros(UInt64, n))

function make_ev(::Type{MyNewMonom}, ev::Vector{U}) where {U} 
    degrees = copy(ev)
    MyNewMonom(degrees)
end

function make_dense!(tmp::Vector{M}, mnm::MyNewMonom) where {M}
    tmp[1:end] = mnm.degrees[1:end]
    tmp
end

function make_hasher(::Type{MyNewMonom}, n::Integer)
    # MonomHash is a global constant type defined in src/types.jl,
    # currently, MonomHash == UInt32
    rand(Groebner.MonomHash, n)
end

function Base.hash(x::MyNewMonom, b::MyNewMonom)
    # plain dot product
    sum(x.degrees .* b.degrees)
end

#=
3. Term order should be defined.
   For now, let's define just one term order, the lexicographic term order.

   - monom_isless(mnm1::MyNewMonom, mnm2::MyNewMonom, ::Lex)
        Returns true if mnm1 < mnm2 in Lex order 
=#

function monom_isless(mnm1::MyNewMonom, mnm2::MyNewMonom, ::Groebner.Lex)
    i = 1
    while i < length(mnm1.degrees) && mnm1.degrees[i] == mnm2.degrees[i]
        i += 1
    end
    return mnm1.degrees[i] < mnm2.degrees[i] ? true : false
end

#=
4. Monomial-by-Monomial arithmetic should be defined.
    This includes:
    - monom_lcm!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
        Writes lcm(a, b) into c, and returns c.
        Assumes c has the same size as a and b.
    - is_gcd_const(a::MyNewMonom, b::MyNewMonom)
        Returns true if gcd(a, b) is constant.
        Assumes the sizes of a and b are the same.
    - monom_product!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
        Writes a * b into c, and returns c.
        Assumes c has the same size as a and b.
    - monom_division!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
        Writes a / b into c, and returns c.
        Assumes c has the same size as a and b.
    - is_monom_divisible(a::MyNewMonom, b::MyNewMonom)
        Returns true if a is divisible by b.
        Assumes the sizes of a and b are the same.
    - is_monom_divisible!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
        Returns true if a is divisible by b, and writes a / b into c.
        Assumes the sizes of a and b are the same.
    - is_monom_elementwise_equal(a::MyNewMonom, b::MyNewMonom)
        Returns true if a and b are elementwise equal.
        Assumes the sizes of a and b are the same.
=#

function monom_lcm!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
    for i in 1:length(c.degrees)
        c.degrees[i] = max(a.degrees[i], b.degrees[i])
    end
    c
end

function is_gcd_const(a::MyNewMonom, b::MyNewMonom)
    for i in 1:length(a.degrees)
        if !iszero(a.degrees[i]) && !iszero(b.degrees[i])
            return false
        end
    end
    return true
end

function monom_product!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
    for i in 1:length(c.degrees)
        c.degrees[i] = a.degrees[i] + b.degrees[i]
    end
    c
end

function monom_division!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
    for i in 1:length(c.degrees)
        c.degrees[i] = a.degrees[i] - b.degrees[i]
    end
    c
end

function is_monom_divisible(a::MyNewMonom, b::MyNewMonom)
    for i in 1:length(a.degrees)
        if a.degrees[i] < b.degrees[i]
            return false
        end
    end
    return true
end

function is_monom_divisible!(c::MyNewMonom, a::MyNewMonom, b::MyNewMonom)
    if is_monom_divisible(a, b)
        monom_division!(c, a, b)
        return true
    end
    return false
end

function is_monom_elementwise_eq(a::MyNewMonom, b::MyNewMonom)
    a.degrees == b.degrees
end

#=
5. A monomial should allow construction of divisibility masks.

   For now, we use a quick hack:
   convert our monomial to a monomial of already existing type PowerVector{UInt64},
   and then compute divisibility mask of it.
=#

function monom_divmask(
        mnm::MyNewMonom,
        DM::Type{Mask},
        ndivvars, divmap,
        ndivbits) where {Mask}
    tmp = Vector{T}(undef, length(mnm.degrees))
    make_dense!(tmp, mnm)
    pv = make_ev(Groebner.PowerVector{UInt64}, tmp)
    Groebner.monom_divmask(pv, DM, ndivvars, divmap, ndivbits)
end
