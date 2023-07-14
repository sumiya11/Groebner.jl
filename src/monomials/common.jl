# Some common utilities for monomial implementations.
#
# This file also provides a description of functions that a monomial type must
# implement

###
# Monomial interface

"""
    construct_const_monom(::Type{Monom}, n)

Constructs a constant monomial of type `Monom` for `n` variables.
"""
function construct_const_monom end

"""
    construct_monom(::Type{Monom}, vector::Vector)

Constructs a monomial of type `Monom` with the partial degrees taken from the
`vector`.
"""
function construct_monom end

"""
    monom_to_dense_vector!(tmp::Vector, monom::Monom)

Returns a dense vector of partial degrees that are stored in `monom`.
Additionaly, writes the result to `tmp`.

"""
function monom_to_dense_vector! end

"""
    totaldeg(monom::Monom)

Returns the total degree of the monomial `monom`.
"""
function totaldeg end

"""
    is_supported_ordering(::Type{Monom}, ord)

Returns `true` if the given monomial type `Monom` implements efficient
comparator functions for the monomial ordering `ord`.
"""
function is_supported_ordering end

"""
    construct_hash_vector(::Type{Monom}, n)

Returns a vector of type Vector{MonomHash}. 

See also the function `monom_hash`.
"""
function construct_hash_vector end

"""
    monom_hash(monom::Monom, hash_vector::Vector)

Retuns a constant of type `MonomHash` that is used for hashing the given
monomial.
"""
function monom_hash end

"""
    max_vars_in_monom(::Type{Monom})

The maximal number of variables that the monomial type `Monom` can potentially
store. 
"""
function max_vars_in_monom end

###
# Common functions

# Threshold for overflow error to trigger
@noinline __monom_overflow_error(c, B) =
    throw(MonomialDegreeOverflow("Overflow may happen with the entry $c of type $B."))

_monom_overflow_threshold(B) = div(typemax(B), 2)

# Checks if there is a risk of exponent overflow
_monom_overflow_check(monom, B) = true
function _monom_overflow_check(deg::Integer, B)
    deg >= _monom_overflow_threshold(B) && __monom_overflow_error(deg, B)
    true
end

# Returns the type of the monomial entry
entrytype(monomtype) = MonomHash
