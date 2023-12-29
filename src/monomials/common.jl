# Some common utilities for monomial implementations.
#
# Provides a description of functions that a monomial type must implement to be
# usable in Groebner. For an example implementation, see the exponent vector in
# monoms/exponentvector.jl
#
# The implementations of the monomial interface must be thread-safe.

###
# Monomial interface

"""
    monom_construct_const_monom(::Type{Monom}, n)

Constructs a constant monomial of type `Monom` capable of storing `n` variables.
"""
function monom_construct_const_monom end

"""
    monom_construct_from_vector(::Type{Monom}, vector::Vector)

Constructs a monomial of type `Monom` with the degrees taken from the given
`vector`.
"""
function monom_construct_from_vector end

"""
    monom_to_vector!(tmp::Vector, monom::Monom)

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

See `monoms/orderings.jl` for details.
"""
function monom_is_supported_ordering end

"""
    monom_construct_hash_vector(::Type{Monom}, n)

Returns a vector of type Vector{MonomHash}. This vector will be later used for
hashing monomials of type `Monom` with `n` variables. 

See also the function `monom_hash`.
"""
function monom_construct_hash_vector end

"""
    monom_hash(monom::Monom, hash_vector::Vector)

Retuns a constant of type `MonomHash` -- the hash of the given `monom`.

It must be that `isequal(monom1, monom2)` implies 
`monom_hash(monom1, hash_vector) == monom_hash(monom2, hash_vector)`. 

The argument `hash_vector` of this function must be constructed using
`monom_construct_hash_vector`.
"""
function monom_hash end

"""
    monom_max_vars(::Type{Monom})

The maximal number of variables that the monomial of type `Monom` can hold. 
"""
function monom_max_vars end

"""
    monom_copy(monom::Monom)

Returns a new monomial of type `Monom` that is equal to the given `monom`.
The returned copy and the original object are independent in memory.
"""
function monom_copy end

###
# Common functions

@noinline __monom_overflow_error(c, B) =
    throw(MonomialDegreeOverflow("Overflow may happen with the entry $c of type $B."))

_monom_overflow_threshold(::Type{T}) where {T <: Integer} = div(typemax(T), 2)

# Checks if there is a risk of exponent overflow
_monom_overflow_check(monom) = true
_monom_overflow_check(monom, B) = true
function _monom_overflow_check(deg::Integer, B)
    deg >= _monom_overflow_threshold(B) && __monom_overflow_error(deg, B)
    true
end

# Returns the type of the monomial entry
monom_entrytype(monomtype) = MonomHash
