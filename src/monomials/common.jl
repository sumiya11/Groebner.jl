# Some common utilities for monomial implementations.
#
# Also provides a description of functions that a monomial type must implement
# to be usable in Groebner. For an example implementation, see the exponent
# vector in monoms/exponentvector.jl

###
# Monomial interface

"""
    construct_const_monom(::Type{Monom}, n)

Constructs a constant monomial of type `Monom` capable of storing `n` variables.
"""
function construct_const_monom end

"""
    construct_monom(::Type{Monom}, vector::Vector)

Constructs a monomial of type `Monom` with the variable degrees read from the
given `vector`.
"""
function construct_monom end

"""
    monom_to_dense_vector!(tmp::Vector, monom::Monom)

Writes the variable degrees of the monomial `monom` to dense vector `tmp` and
returns `tmp`.
"""
function monom_to_dense_vector! end

"""
    totaldeg(monom::Monom)

Returns the total degree of the monomial `monom`.
"""
function totaldeg end

"""
    is_supported_ordering(::Type{Monom}, ord<:AbstractMonomialOrdering)

Returns `true` if the given monomial type `Monom` implements efficient
comparator function for the monomial ordering `ord`.

See `monoms/orderings.jl` for details.
"""
function is_supported_ordering end

"""
    construct_hash_vector(::Type{Monom}, n)

Returns a dense vector of type Vector{MonomHash}. This vector will be later used
for hashing monomials of type `Monom` with `n` variables. 

See also the function `monom_hash`.
"""
function construct_hash_vector end

"""
    monom_hash(monom::Monom, hash_vector::Vector)

Retuns a constant of type `MonomHash` --- the hash of the given `monom`.

The argument `hash_vector` will be constructed using the function
`construct_hash_vector`.
"""
function monom_hash end

"""
    max_vars_in_monom(::Type{Monom})

The maximal number of variables that the monomial type `Monom` can hold. 
"""
function max_vars_in_monom end

"""
    copy_monom(monom)

Returns a new monom that is equal to `monom`.
"""
function copy_monom end

###
# Common functions

# Threshold for overflow error to trigger
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
