
function column_to_monom_mapping!(graph, matrix, symbol_ht)
    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        if hdata[i].idx == 2
            k += 1
        end
    end

    column_to_monom = matrix.column_to_monom

    matrix.nleft = k  # CHECK!
    # -1 as long as hashtable load is always 1 more than actual
    matrix.nright = load - matrix.nleft - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        hdata[column_to_monom[k]].idx = k
    end

    @inbounds for k in 1:(matrix.nupper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nlower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    matrix.ncolumns = matrix.nleft + matrix.nright

    @assert matrix.nleft + matrix.nright == symbol_ht.load - 1 == matrix.ncolumns
    @assert matrix.nlower + matrix.nupper == matrix.nrows
end

function reduction_apply!(
    graph::ComputationGraphF4,
    ring,
    basis,
    matrix,
    ht,
    rng,
    symbol_ht,
    iter
)
    if length(graph.matrix_sorted_columns) < iter
        column_to_monom_mapping!(matrix, symbol_ht)
        push!(graph.matrix_sorted_columns, matrix.column_to_monom)
    else
        matrix.column_to_monom = graph.matrix_sorted_columns[iter]
        column_to_monom_mapping!(graph, matrix, symbol_ht)
    end

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    @log level = -3 repr_matrix(matrix)

    @log level = -6 "Apply: after mapping and sorting" matrix
    flag = linear_algebra!(graph, ring, matrix, basis, :apply, rng)
    if !flag
        return false
    end
    @log level = -6 "After linear algebra" matrix

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
    true
end

function symbolic_preprocessing!(
    graph::ComputationGraphF4,
    iters::Int,
    basis,
    matrix,
    hashtable,
    symbol_ht
)
    lowrows, lowmults = graph.matrix_lower_rows[iters]
    uprows, upmults = graph.matrix_upper_rows[iters]
    matrix_info = graph.matrix_infos[iters]
    nonzeroed_rows = graph.matrix_nonzeroed_rows[iters]
    nlow = length(nonzeroed_rows) # matrix_info.nlow # length(nonzeroed_rows)
    nup = length(uprows) # matrix_info.nup # length(uprows)

    matrix.upper_rows = Vector{Vector{ColumnLabel}}(undef, nup)
    matrix.lower_rows = Vector{Vector{ColumnLabel}}(undef, nlow)
    matrix.lower_to_coeffs = Vector{Int}(undef, nlow)
    matrix.upper_to_coeffs = Vector{Int}(undef, nup)

    resize_hashtable_if_needed!(symbol_ht, nlow + nup + 2)
    for i in 1:nlow
        mult_idx = lowmults[i]
        poly_idx = lowrows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        # vidx = insert_in_hash_table!(symbol_ht, etmp)

        matrix.lower_rows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)
        symbol_ht.hashdata[matrix.lower_rows[i][1]].idx = PIVOT_COLUMN

        matrix.lower_to_coeffs[i] = poly_idx
    end

    for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]
        # vidx = insert_in_hash_table!(symbol_ht, etmp)

        matrix.upper_rows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)
        symbol_ht.hashdata[matrix.upper_rows[i][1]].idx = PIVOT_COLUMN

        matrix.upper_to_coeffs[i] = poly_idx
    end

    i = MonomIdx(symbol_ht.offset)
    @inbounds while i <= symbol_ht.load
        # not a reducer
        if iszero(symbol_ht.hashdata[i].idx)
            symbol_ht.hashdata[i].idx = UNKNOWN_PIVOT_COLUMN
        end
        i += MonomIdx(1)
    end

    @log level = -6 "Symbol ht:" symbol_ht

    matrix.nrows = nlow + nup
    matrix.nlower = nlow
    matrix.nupper = nup
    matrix.size = matrix.nrows
end

function reducegb_f4_apply!(
    graph,
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    hashtable::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    iter
) where {M}
    @log level = -5 "Entering apply autoreduction" basis

    lowrows, lowmults = graph.matrix_lower_rows[end]
    uprows, upmults = graph.matrix_upper_rows[end]
    matrix_info = graph.matrix_infos[end]
    nonzeroed_rows = graph.matrix_nonzeroed_rows[end]
    nlow = length(lowrows) # matrix_info.nlow
    nup = length(uprows) # matrix_info.nup # 
    # ncols = matrix_info.ncols
    nonred = graph.nonredundant_indices_before_reduce

    matrix.upper_rows = Vector{Vector{ColumnLabel}}(undef, nup)
    matrix.lower_rows = Vector{Vector{ColumnLabel}}(undef, nlow)
    matrix.lower_to_coeffs = Vector{Int}(undef, nlow)
    matrix.upper_to_coeffs = Vector{Int}(undef, nup)

    matrix.nrows = 0
    matrix.ncolumns = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nupper = nup
    matrix.nlower = nlow

    etmp = construct_const_monom(M, hashtable.nvars)
    # etmp is now set to zero, and has zero hash

    resize_hashtable_if_needed!(symbol_ht, nlow + nup + 2)

    # needed for correct column count in symbol hashtable
    matrix.ncolumns = matrix.nrows
    matrix.nupper = matrix.nrows

    @log level = -5 "Before autoreduce apply" basis uprows upmults lowmults matrix_info

    for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        @log level = -6 "Up row" etmp

        # Not too good
        # iszero(h) && continue

        # vidx = insert_in_hash_table!(symbol_ht, etmp)
        matrix.upper_rows[i] =
            multiplied_poly_to_matrix_row!(symbol_ht, hashtable, h, etmp, rpoly)

        matrix.upper_to_coeffs[i] = poly_idx
    end

    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.hashdata[i].idx = UNKNOWN_PIVOT_COLUMN
    end
    # matrix.ncolumns = ncols

    matrix.nrows = nlow + nup
    matrix.nlower = nlow
    matrix.nupper = nup
    matrix.size = matrix.nrows

    @log level = -2 length(graph.matrix_sorted_columns) iter
    if length(graph.matrix_sorted_columns) < iter
        column_to_monom_mapping!(matrix, symbol_ht)
        push!(graph.matrix_sorted_columns, matrix.column_to_monom)
    else
        matrix.column_to_monom = graph.matrix_sorted_columns[iter]
        column_to_monom_mapping!(graph, matrix, symbol_ht)
    end
    matrix.ncolumns = matrix.nleft + matrix.nright

    sort_matrix_upper_rows_decreasing!(matrix)

    @log level = -6 "In autoreduction apply" basis matrix

    flag = deterministic_sparse_rref_interreduce!(ring, matrix, basis)
    if !flag
        return false
    end
    convert_rows_to_basis_elements!(matrix, basis, hashtable, symbol_ht)

    basis.nfilled = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    @log level = -6 "Apply autoreduction: after linalg" basis matrix

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    output_nonredundant = graph.output_nonredundant_indices
    # basis.nonredundant = output_nonredundant
    for i in 1:length(output_nonredundant)
        basis.nonredundant[i] = output_nonredundant[i]
        basis.divmasks[i] =
            hashtable.hashdata[basis.monoms[basis.nonredundant[i]][1]].divmask
    end
    basis.nnonredundant = length(output_nonredundant)
    true
end

function standardize_basis_apply!(graph)
    basis = graph.gb_basis
    buf = graph.buf_basis
    basis.size = basis.nprocessed = basis.nfilled = basis.nnonredundant
    output_nonredundant = graph.output_nonredundant_indices
    output_sort = graph.output_sort_indices
    @inbounds for i in 1:(basis.nprocessed)
        basis.coeffs[i] = buf.coeffs[output_nonredundant[output_sort[i]]]
    end
    buf.nprocessed = buf.nnonredundant = 0
    buf.nfilled = graph.input_basis.nfilled
    normalize_basis!(graph.ring, basis)
end

function standardize_basis_learn!(graph, ring, basis, ht, ord)
    @inbounds for i in 1:(basis.nnonredundant)
        idx = basis.nonredundant[i]
        basis.nonredundant[i] = i
        basis.isredundant[i] = false
        basis.coeffs[i] = basis.coeffs[idx]
        basis.monoms[i] = basis.monoms[idx]
    end
    basis.size = basis.nprocessed = basis.nfilled = basis.nnonredundant
    resize!(basis.coeffs, basis.nprocessed)
    resize!(basis.monoms, basis.nprocessed)
    resize!(basis.divmasks, basis.nprocessed)
    resize!(basis.nonredundant, basis.nprocessed)
    resize!(basis.isredundant, basis.nprocessed)
    perm = sort_polys_by_lead_increasing!(basis, ht, ord=ord)
    graph.output_sort_indices = perm
    normalize_basis!(ring, basis)
end

function f4_apply!(graph, ring, basis::Basis{C}, params) where {C <: Coeff}
    @invariant basis_well_formed(:input_f4_apply!, ring, basis, graph.hashtable)
    @assert params.reduced == true

    normalize_basis!(ring, basis)

    iters_total = length(graph.matrix_infos) - 1
    iters = 0
    hashtable = graph.hashtable

    matrix = initialize_matrix(ring, C)
    @log level = -6 "Applying modulo $(ring.ch)"

    update_basis!(basis, hashtable)

    @log level = -6 "Input basis:" basis

    while iters < iters_total
        iters += 1
        @log level = -3 "F4 Apply iteration $iters"

        symbol_ht = initialize_secondary_hashtable(hashtable)

        symbolic_preprocessing!(graph, iters, basis, matrix, hashtable, symbol_ht)
        @log level = -5 "After symbolic preprocessing:" matrix

        flag = reduction_apply!(
            graph,
            ring,
            basis,
            matrix,
            hashtable,
            params.rng,
            symbol_ht,
            iters
        )
        if !flag
            # Unlucky cancellation of basis coefficients happened
            return false
        end

        @log level = -6 "After reduction_apply:" matrix basis

        update_basis!(basis, hashtable)

        @show_locals

        @log level = -6 "After update apply" basis

        matrix = initialize_matrix(ring, C)
    end

    @log level = -6 "Before reduction" basis

    # mark_redundant!(basis)

    if params.reduced
        @log level = -6 "Autoreducing the final basis.."
        symbol_ht = initialize_secondary_hashtable(hashtable)
        flag = reducegb_f4_apply!(
            graph,
            ring,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            iters_total + 1
        )
        if !flag
            return false
        end
        @log level = -6 "Autoreduced!"
    end

    standardize_basis_apply!(graph)
    basis = graph.gb_basis

    @log level = -6 "After apply standardization" basis

    @invariant basis_well_formed(:output_f4_apply!, ring, basis, hashtable)

    true
end

function update!(graph, pairset, basis, hashtable, update_ht)
    update!(pairset, basis, hashtable, update_ht)
end

function symbolic_preprocessing!(graph, basis, matrix, hashtable, symbol_ht)
    symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
    @log level = -6 "Symbol ht:" symbol_ht
end

function reduction_learn!(graph, ring, basis, matrix, hashtable, symbol_ht, linalg, rng)
    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    @log level = -6 "Learn: after mapping and sorting" matrix

    @log level = -3 repr_matrix(matrix)

    linear_algebra!(graph, ring, matrix, basis, :learn, rng)

    @log level = -6 "After linear algebra" matrix

    convert_rows_to_basis_elements!(matrix, basis, hashtable, symbol_ht)
end

function f4_reducegb_learn!(
    graph,
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M}
) where {M}
    @log level = -6 "Entering learn autoreduction" basis

    #=

    =#

    etmp = construct_const_monom(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    reinitialize_matrix!(matrix, basis.nnonredundant)
    uprows = matrix.upper_rows

    @log level = -6 "Before autoreduce learn" basis

    # add all non redundant elements from basis
    # as matrix upper rows
    @inbounds for i in 1:(basis.nnonredundant) #
        matrix.nrows += 1
        uprows[matrix.nrows] = multiplied_poly_to_matrix_row!(
            symbol_ht,
            ht,
            MonomHash(0),
            etmp,
            basis.monoms[basis.nonredundant[i]]
        )

        matrix.upper_to_coeffs[matrix.nrows] = basis.nonredundant[i]
        matrix.upper_to_mult[matrix.nrows] = insert_in_hash_table!(ht, etmp)
        # set lead index as 1
        symbol_ht.hashdata[uprows[matrix.nrows][1]].idx = UNKNOWN_PIVOT_COLUMN
    end
    graph.nonredundant_indices_before_reduce = basis.nonredundant[1:(basis.nnonredundant)]

    # needed for correct column count in symbol hashtable
    matrix.ncolumns = matrix.nrows
    matrix.nupper = matrix.nrows

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.hashdata[i].idx = UNKNOWN_PIVOT_COLUMN
    end

    column_to_monom_mapping!(matrix, symbol_ht)
    matrix.ncolumns = matrix.nleft + matrix.nright

    sort_matrix_upper_rows_decreasing!(matrix)

    @log level = -6 "In autoreduction learn" basis matrix

    deterministic_sparse_rref_interreduce_learn!(graph, ring, matrix, basis)

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)

    basis.nfilled = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    @log level = -6 "Learn autoreduction: after linalg" basis matrix

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    k = 0
    i = 1
    @label Letsgo
    @inbounds while i <= basis.nprocessed
        @inbounds for j in 1:k
            if is_monom_divisible(
                basis.monoms[basis.nfilled - i + 1][1],
                basis.monoms[basis.nonredundant[j]][1],
                ht
            )
                i += 1
                @goto Letsgo
            end
        end
        k += 1
        basis.nonredundant[k] = basis.nfilled - i + 1
        basis.divmasks[k] = ht.hashdata[basis.monoms[basis.nonredundant[k]][1]].divmask
        i += 1
    end
    basis.nnonredundant = k

    graph.output_nonredundant_indices = copy(basis.nonredundant[1:k])
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

    @assert basis == graph.gb_basis
    @assert params.reduced === true

    @log level = -3 "Entering F4 Learn phase."
    normalize_basis!(ring, basis)

    matrix = initialize_matrix(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = initialize_secondary_hashtable(hashtable)
    symbol_ht = initialize_secondary_hashtable(hashtable)

    # add the first batch of critical pairs to the pairset
    @log level = -3 "Processing initial polynomials, generating first critical pairs"
    pairset_size = update!(graph, pairset, basis, hashtable, update_ht)
    @log level = -3 "Out of $(basis.nfilled) polynomials, $(basis.nprocessed) are non-redundant"
    @log level = -3 "Generated $(pairset.load) critical pairs"

    @log level = -6 "Input basis:" basis

    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log level = -3 "F4: iteration $i"
        @log level = -3 "F4: available $(pairset.load) pairs"

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
        @log level = -6 "Formed a matrix of size X, DISPLAY_MATRIX"
        @log level = -6 "After symbolic preprocessing:" matrix

        # reduces polys and obtains new potential basis elements
        reduction_learn!(
            graph,
            ring,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            params.linalg,
            params.rng
        )

        @log level = -6 "After reduction_learn:" matrix basis

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = update!(graph, pairset, basis, hashtable, update_ht)

        @log level = -6 "After update learn" basis
        # clear symbolic hashtable
        # clear matrix
        matrix    = initialize_matrix(ring, C)
        symbol_ht = initialize_secondary_hashtable(hashtable)

        if i > 10_000
            @log level = 1 "Something has gone wrong in F4. Error will follow."
            @show_locals
            __throw_maximum_iterations_exceeded(i)
        end
    end

    @log level = -6 "Before filter redundant" basis

    if params.sweep
        @log level = -3 "Sweeping redundant elements in the basis"
        sweep_redundant!(basis, hashtable)
    end

    # mark redundant elements
    mark_redundant!(basis)
    @log level = -3 "Filtered elements marked redundant"

    @log level = -6 "Before autoreduction" basis

    if params.reduced
        @log level = -2 "Autoreducing the final basis.."
        f4_reducegb_learn!(graph, ring, basis, matrix, hashtable, symbol_ht)
        @log level = -3 "Autoreduced!"
    end

    @log level = -3 "Finalizing computation graph"
    finalize_graph!(graph)

    standardize_basis_learn!(graph, ring, basis, hashtable, hashtable.ord)

    @log level = -6 "After learn standardization" basis

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    @invariant basis_well_formed(:output_f4_learn!, ring, basis, hashtable)

    nothing
end
