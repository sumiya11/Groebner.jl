# Input-output conversions of polynomials.
# Currently, conversions work with polynomials from
#  - AbstractAlgebra.jl
#  - DynamicPolynomials.jl
#  - Nemo.jl
#  - Singular.jl (TODO: add tests)

# TODO: move to input-output.jl !
    # iszerox = remove_zeros_from_input!(ring, exps, coeffs)
    # iszerox && (return convert_to_output(ring, polynomials, exps, coeffs, metainfo))
    
#=
    Our conventions:
    - Trying to compute a Groebner basis of an empty set is an error.
    - The Groebner basis of [0] is [0].
    - The Groebner basis of [0,..., 0] is [0].
    - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]
=#

#=
    A note on how we represent polynomials.

    First, coefficients, exponents and polynomial ring information
    are extracted from input polynomials with `convert_to_internal`.

    Inside the algorithm all monomials are hashed without collisions,
    so that an integer represents a single monomial.
    A single monomial could be stored in several different ways, 
    according to the prescibed monomial implementation.

    A polynomial is represented with a dynamic coefficients vector
    together with a dynamic vector of hashtable indices of monomials.

    After the basis is computed, hash table indices are converted back to 
    monomials, and the `convert_to_output` is called 
    to convert internal structures to the original polynomial type.
=#

#=
    A note about monomial orderings.

    Groebner.jl supports several monomial orderings 
    (all of which are subtypes of AbstractMonomialOrdering, 
    see src/monoms/orderings.jl for details)

    Polynomials from AbstractAlgebra.jl, DynamicPolynomials.jl, 
    and some other packages, do not support some of the orderings supported by Groebner.jl.

    We compute the basis in the requested ordering in Groebner.jl,
    and then output the polynomials in the ordering of the input.

    For example, say that the input polynomials are from AbstractAlgebra.jl 
    and are in the Lex ordering. The requested basis is in some Weighted ordering.
    Then, the output should be a correct basis in the Weighted ordering,
    though the terms in the output are ordered w.r.t. Lex ordering.
=#

# For compatibility with polynomials from AbstractAlgebra.jl and Nemo.jl
const _AA_supported_orderings_symbols = (:lex, :deglex, :degrevlex)
# The type of exponent vector entries used internally in AbstractAlgebra.jl
const _AA_exponenttype = UInt64

"""
    Contains information about some polynomial ring.
"""
mutable struct PolyRing{Char <: CoeffFF, Ord <: AbstractMonomialOrdering}
    # number of variables
    nvars::Int
    # monomial ordering
    ord::Ord
    # characteristic of the coefficient field
    ch::Char
    # information about the original ring of input. Options are:
    #    :AbstractAlgebra for AbstractAlgebra,
    #    :MultivariatePolynomials for subtypes of MultivariatePolynomials.jl, 
    #       e.g, for DynamicPolynomials.jl,
    #    :hasparent for polynomials constructed with parent ring, e.g., Nemo
    #    :undefined
    origring::Symbol
end

"""
    Converts input polynomials to internal representation used by the algorithm.
    Extracts polynomial ring information, and polynomial exponents and coefficients.

    This is the most general implementation.
    Works for polynomials that implement AbstractAlgebra.Generic.MPoly intefrace:
        . `AbstractAlgebra`
        . `Nemo`
        . `Singular`

    Currently, there are also specializations for polynomials from:
        . `MultivariatePolynomials`

"""
function convert_to_internal(
    representation,
    orig_polys::Vector{T},
    ordering::AbstractMonomialOrdering
) where {T}
    isempty(orig_polys) && throw(DomainError(orig_polys, "Empty input."))

    if hasmethod(AbstractAlgebra.parent, Tuple{typeof(first(orig_polys))})
        convert_to_internal(representation, orig_polys, ordering, Val(:hasparent))
    else
        convert_to_internal(representation, orig_polys, ordering, Val(:undefined))
    end
end

function convert_to_internal(
    representation,
    orig_polys::Vector{T},
    ordering::AbstractMonomialOrdering,
    ::Val{:undefined}
) where {T}
    error(
        "Sorry, we don't work with this type of polynomials ($T) yet. Feel free to open an issue"
    )
end
