# This file is a part of Groebner.jl. License is GNU GPL v2.

###
#  Autoreducing polynomials

function autoreduce2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return [monoms[1]], [coeffs[1]]
    end
    monoms, coeffs = _monoms, _coeffs
    _autoreduce2(ring, monoms, coeffs, params)
end

function _autoreduce2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    basis, _, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    basis_update!(basis, hashtable)
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(hashtable)
    f4_autoreduce!(ring, basis, matrix, hashtable, symbol_ht, params)
    basis_standardize!(ring, basis, hashtable, params.arithmetic, params.changematrix)
    basis_export_data(basis, hashtable)
end
