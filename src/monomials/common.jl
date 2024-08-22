# This file is a part of Groebner.jl. License is GNU GPL v2.

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
