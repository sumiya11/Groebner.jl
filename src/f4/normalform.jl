

function normal_form_f4(
            basis::Basis,
            ht::MonomialHashtable,
            tobereduced::Basis)

    matrix = initialize_matrix()
    symbol_ht = initialize_secondary_hash_table(ht)

    select_tobereduced!(tobereduced, matrix, symbol_ht, ht)
    symbolic_preprocessing!(matrix, basis, symbol_ht, ht)

    convert_hashes_to_columns!(matrix, symbol_ht)

    sort_matrix_rows_decreasing!(matrix) # for pivots,  AB part

    exact_sparse_linear_algebra_nf!(matrix, tobereduced, basis)

    convert_nf_rows_to_basis_elements!(matrix, tobereduced, ht, symbol_ht)

    tobereduced
end
