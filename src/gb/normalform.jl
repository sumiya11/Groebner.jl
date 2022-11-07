
#------------------------------------------------------------------------------

function clean_input_normalform!(ring::PolyRing, 
        basisexps, basiscoeffs, tbrexps, tbrcoeffs)
    clean_input_groebner!(ring, basisexps, basiscoeffs)
end

#------------------------------------------------------------------------------

function normal_form_f4!(
            ring::PolyRing,
            basis::Basis{C},
            ht::MonomialHashtable,
            tobereduced::Basis{C}) where {C<:Coeff}

    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hash_table(ht)

    select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)

    convert_hashes_to_columns!(matrix, symbol_ht)

    #dump(matrix, maxdepth=5)
    #@error "((()))"

    sort_matrix_rows_decreasing!(matrix) # for pivots,  AB part

    #=
    dump(matrix, maxdepth=5)
    @error "UPROWS"
    for (i, u) in enumerate(matrix.uprows)
        print(map(c -> symbol_ht.exponents[matrix.col2hash[c]], u))
        println(" ", basis.coeffs[matrix.up2coef[i]])
    end
    @error "LOWrows"
    for (i,u) in enumerate(matrix.lowrows)
        println(map(c -> symbol_ht.exponents[matrix.col2hash[c]], u))
        println(" ", tobereduced.coeffs[matrix.low2coef[i]])
    end

    @error "" Int(ring.ch)
    =#

    exact_sparse_linear_algebra_nf!(ring, matrix, tobereduced, basis)

    # @warn "lll"
    # dump(matrix, maxdepth=5)

    convert_nf_rows_to_basis_elements!(matrix, tobereduced, ht, symbol_ht)

    tobereduced
end

function normal_form_f4(
                ring::PolyRing,
                basisexps::Vector{Vector{M}},
                basiscoeffs::Vector{Vector{C}},
                tobereducedexps::Vector{Vector{M}},
                tobereducedcfs::Vector{Vector{C}},
                rng::Rng,
                tablesize::Int=2^16) where {M, Rng, C<:Coeff}

    basis, ht = initialize_structures(ring, basisexps,
                                        basiscoeffs, rng, tablesize)

    tobereduced, ht = initialize_structures_nf(ring, tobereducedexps,
                                        tobereducedcfs, rng, tablesize, ht)
                                        
    tobereduced = normal_form_f4!(ring, basis, ht, tobereduced)

    export_basis_data(tobereduced, ht)
end
