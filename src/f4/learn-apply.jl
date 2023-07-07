
function reduction!(graph::ComputationGraphF4, ring, basis, matrix, ht)
    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    linear_algebra!(graph, ring, matrix, basis, linalg, rng)

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
end

function symbolic_preprocessing!(graph::ComputationGraphF4, basis, matrix, hashtable) end

function f4_apply!(graph, ring, basis::Basis{C}, params) where {M <: Monom, C <: Coeff}
    @invariant basis_well_formed(:input_f4_apply!, ring, basis, hashtable)

    normalize_basis!(ring, basis)

    iters_total = length(graph.matrix_infos)
    iters = 0
    hashtable = graph.hashtable

    matrix = initialize_matrix(ring, C)

    while iters < iters_total
        @log "F4 Apply iteration $iters"

        symbolic_preprocessing!(graph, basis, matrix, hashtable)
        reduction!(graph, ring, basis, matrix, hashtable)

        matrix = initialize_matrix(ring, C)
        iters += 1
    end

    if params.reduced
        @log "Autoreducing the final basis.."
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht)
        @log "Autoreduced!"
    end
    nothing
end

function update!(graph, pairset, basis, hashtable, update_ht)
    update!(pairset, basis, hashtable, update_ht)
end

function symbolic_preprocessing!(graph, basis, matrix, hashtable, symbol_ht)
    symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
end

function reduction!(graph, ring, basis, matrix, hashtable, symbol_ht, linalg, rng)
    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    linear_algebra!(graph, ring, matrix, basis, linalg, rng)

    convert_rows_to_basis_elements!(matrix, basis, hashtable, symbol_ht)
end

function f4_learn!(
    graph,
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    # @invariant hashtable_well_formed(:input_f4!, ring, hashtable)
    @invariant basis_well_formed(:input_f4_learn!, ring, basis, hashtable)
    # @invariant pairset_well_formed(:input_f4!, pairset, basis, ht)

    @log level = 5 "Entering F4 Learn phase."
    # TODO: decide on the number field arithmetic implementation
    normalize_basis!(ring, basis)

    matrix = initialize_matrix(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = initialize_secondary_hashtable(hashtable)
    symbol_ht = initialize_secondary_hashtable(hashtable)

    # add the first batch of critical pairs to the pairset
    @log level = 6 "Processing initial polynomials, generating first critical pairs"
    pairset_size = update!(graph, pairset, basis, hashtable, update_ht)
    @log level = 6 "Out of $(basis.ntotal) polynomials, $(basis.nprocessed) are non-redundant"
    @log level = 6 "Generated $(pairset.load) critical pairs"

    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log "F4: iteration $i"
        @log "F4: available $(pairset.load) pairs"

        # selects pairs for reduction from pairset following normal strategy
        # (minimal lcm degrees are selected),
        # and puts these into the matrix rows
        select_normal!(
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            maxpairs=params.maxpairs
        )
        # Color with [F4]

        symbolic_preprocessing!(graph, basis, matrix, hashtable, symbol_ht)
        @log "Formed a matrix of size X, DISPLAY_MATRIX"

        # reduces polys and obtains new potential basis elements
        reduction!(
            graph,
            ring,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            params.linalg,
            params.rng
        )
        # TODO: insert monoms from symbol_ht into the main hashtable
        @log ""

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = update!(graph, pairset, basis, hashtable, update_ht)
        # TODO: move the above
        @log "something something"

        # clear symbolic hashtable
        # clear matrix
        matrix    = initialize_matrix(ring, C)
        symbol_ht = initialize_secondary_hashtable(hashtable)

        if i > 10_000
            # TODO: log useful info here
            @log level = 100 "Something has probably gone wrong in F4. An exception will follow."
            __error_maximal_number_exceeded(
                "Something has probably gone wrong in F4. Please submit a github issue."
            )
        end
    end

    # remove redundant elements
    # filter_redundant!(basis)
    # @log "Filtered elements marked redundant"

    if params.reduced
        @log "Autoreducing the final basis.."
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht)
        @log "Autoreduced!"
    end

    # standardize_basis!(ring, basis, hashtable, hashtable.ord)

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    # @invariant basis_well_formed(:output_f4!, ring, basis, hashtable)

    nothing
end
