
# autoreduce for Finite fields
function _autoreduce(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    basis, _, hashtable = initialize_structs(ring, monoms, coeffs, params)
    update_basis!(basis, hashtable)
    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hashtable(hashtable)
    reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht, params)
    standardize_basis!(ring, basis, hashtable, hashtable.ord, params.arithmetic)
    export_basis_data(basis, hashtable)
end
