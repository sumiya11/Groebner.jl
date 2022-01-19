
function select_tobereduced!(
                basis::Basis, tobereduced::Basis,
                matrix::MacaulayMatrix,
                symbol_ht::MonomialHashtable, ht::MonomialHashtable)

    # prepare to load all elems from tobereduced
    # into low rows
    reinitialize_matrix!(matrix, 2^6)
    resize!(matrix.lowrows, tobereduced.ntotal)

    # TODO
    etmp = zeros(UInt16, ht.explen)

    for i in 1:tobereduced.ntotal
        matrix.nrows += 1
        gen = tobereduced.gens[i]
        h = UInt32(0)
        matrix.lowrows[matrix.nrows] = multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, gen)
        matrix.low2coef[matrix.nrows] = i
    end

    basis.ntotal
    basis.nlead = basis.ndone = basis.ntotal
    basis.isred .= 0
    for i in 1:basis.nlead
        basis.nonred[i] = i
        basis.lead[i] = ht.hashdata[basis.gens[i][1]].divmask
    end
end

function normal_form_f4!(
            ring::PolyRing,
            basis::Basis,
            ht::MonomialHashtable,
            tobereduced::Basis)

    matrix = initialize_matrix(ring)
    symbol_ht = initialize_secondary_hash_table(ht)

    select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)

    convert_hashes_to_columns!(matrix, symbol_ht)

    sort_matrix_rows_decreasing!(matrix) # for pivots,  AB part
    exact_sparse_linear_algebra_nf!(matrix, tobereduced, basis)

    convert_nf_rows_to_basis_elements!(matrix, tobereduced, ht, symbol_ht)

    tobereduced
end

function normal_form_f4(
                ring::PolyRing,
                basisexps::Vector{Vector{Vector{UInt16}}},
                basiscoeffs::Vector{Vector{UInt64}},
                tobereducedexps::Vector{Vector{Vector{UInt16}}},
                tobereducedcfs::Vector{Vector{UInt64}},
                rng::Rng,
                tablesize::Int=2^16) where {Rng}

    basis, ht = initialize_structures(ring, basisexps,
                                        basiscoeffs, rng, tablesize)

    tobereduced, ht = initialize_structures(ring, tobereducedexps,
                                        tobereducedcfs, rng, tablesize,
                                        ht)
    tobereduced = normal_form_f4!(ring, basis, ht, tobereduced)

    export_basis_data(tobereduced, ht)
end
