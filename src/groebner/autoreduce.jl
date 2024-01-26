# This file is a part of Groebner.jl. License is GNU GPL v2.

###
#  Autoreducing a set of polynomials

function _autoreduce1(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    basis, _, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    basis_update!(basis, hashtable)
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(hashtable)
    f4_autoreduce!(ring, basis, matrix, hashtable, symbol_ht, params)
    basis_standardize!(ring, basis, hashtable, hashtable.ord, params.arithmetic)
    basis_export_data(basis, hashtable)
end
