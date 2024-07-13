# This file is a part of Groebner.jl. License is GNU GPL v2.

# Common utilities for monomial implementations. 
#
# Also describes the monomial interface, that is, the list of functions that a
# monomial type must implement in order to be usable in Groebner.

###
# Monomial interface

"""
    monom_construct_const(::Type{Monom}, n)

Constructs a constant monomial of type `Monom` capable of storing `n` variables.
"""
function monom_construct_const end

"""
    monom_construct_from_vector(::Type{Monom}, vector::AbstractVector)

Constructs a monomial of type `Monom` with the degrees taken from the given
`vector`.
"""
function monom_construct_from_vector end

"""
    monom_to_vector!(tmp::AbstractVector, monom::Monom)

Writes the degrees of the monomial `monom` to vector `tmp` and returns `tmp`.
"""
function monom_to_vector! end

"""
    monom_totaldeg(monom::Monom)

Returns the total degree of the monomial `monom`.
"""
function monom_totaldeg end

"""
    monom_is_supported_ordering(::Type{Monom}, ord <: AbstractMonomialOrdering)

Returns `true` if the given monomial type `Monom` implements efficient
comparator function for the monomial ordering `ord`.
"""
function monom_is_supported_ordering end

"""
    monom_construct_hash_vector(rng, ::Type{Monom}, n)

Returns a vector of type Vector{MonomHash}. This vector will be later used for
hashing monomials of type `Monom` with `n` variables via `monom_hash`.
"""
function monom_construct_hash_vector end

"""
    monom_hash(monom::Monom, hash_vector::Vector)

Retuns the hash of the given `monom` using the given `hash_vector`.

`isequal(monom1, monom2)` must imply
`monom_hash(monom1, hash_vector) == monom_hash(monom2, hash_vector)`. 

The argument `hash_vector` of this function must be constructed using
`monom_construct_hash_vector`.
"""
function monom_hash end

"""
    monom_max_vars(::Type{Monom})

The number of variables that the monomial of type `Monom` can possibly hold. 
"""
function monom_max_vars end

"""
    monom_copy(monom::Monom)

Returns a copy of `monom`.
The returned copy and the original object must not alias in memory.
"""
function monom_copy end

###
# Common functions

"""
    _monom_overflow_check(monom::Monom)

Checks if there is a risk of monomial exponent overflow. 
    
If overflow if possible, throws a `MonomialDegreeOverflow`.
"""
function _monom_overflow_check end

@noinline __monom_overflow_error(c, B) =
    throw(MonomialDegreeOverflow("Overflow may happen with the entry $c of type $B."))

_monom_overflow_threshold(::Type{T}) where {T <: Integer} = div(typemax(T), 2)

_monom_overflow_check(monom) = true
_monom_overflow_check(monom, B) = true
function _monom_overflow_check(deg::Integer, B)
    deg >= _monom_overflow_threshold(B) && __monom_overflow_error(deg, B)
    true
end

# Returns the type of the monomial entry
monom_entrytype(monom) = MonomHash
