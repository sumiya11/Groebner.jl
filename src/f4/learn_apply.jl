# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# F4 learn & apply with tracing
#
# Tracing is implemented as an add-on to the F4 algorithm. Here, we try to reuse
# the functions from F4 as much as possible. At the same time, F4 is completely
# independent from this file.

###
# F4 learn stage

function f4_initialize_structs_with_trace(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    make_monic=true,
    sort_input=true
) where {M <: Monom, C <: Coeff}
    basis, pairset, hashtable, permutation = f4_initialize_structs(
        ring,
        monoms,
        coeffs,
        params,
        make_monic=make_monic,
        sort_input=sort_input
    )

    trace = trace_initialize(ring, basis_deepcopy(basis), basis, hashtable, permutation, params)

    trace, basis, pairset, hashtable, permutation
end

function matrix_compute_pivot_signature(pivots::Vector{Vector{MonomId}}, from::Int, sz::Int)
    sgn = UInt64(0x7e2d6fb6448beb77)
    sgn = sgn - UInt64(89 * sz)
    @inbounds for i in from:(from + sz - 1)
        p = pivots[i]
        sgni = zero(UInt64)
        for j in 1:length(p)
            sgni = p[j] - UInt64(13) * sgni
        end
        sgn = sgn - UInt(13) * sgni
    end
    sgn
end

function f4_reduction_learn!(
    trace::Trace,
    basis::Basis,
    matrix::MacaulayMatrix,
    hashtable::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters
)
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    linalg_main!(matrix, basis, params, trace, linalg=LinearAlgebra(:learn, :sparse))
    matrix_convert_rows_to_basis_elements!(
        matrix,
        basis,
        hashtable,
        symbol_ht,
        params;
        batched_ht_insert=true
    )
    pivot_indices = map(i -> Int32(basis.monoms[basis.n_processed + i][1]), 1:(matrix.npivots))
    push!(trace.matrix_pivot_indices, pivot_indices)
    matrix_pivot_signature =
        matrix_compute_pivot_signature(basis.monoms, basis.n_processed + 1, matrix.npivots)
    push!(trace.matrix_pivot_signatures, matrix_pivot_signature)
end

function f4_reducegb_learn!(
    trace::Trace,
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom}
    etmp = monom_construct_const(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    matrix_reinitialize!(matrix, basis.n_nonredundant)
    uprows = matrix.upper_rows

    # add all non redundant elements from basis
    # as matrix upper rows
    @inbounds for i in 1:(basis.n_nonredundant) #
        row_idx = matrix.nrows_filled_upper += 1
        uprows[row_idx] = matrix_polynomial_multiple_to_row!(
            matrix,
            symbol_ht,
            ht,
            MonomHash(0),
            etmp,
            basis.monoms[basis.nonredundant_indices[i]]
        )

        matrix.upper_to_coeffs[row_idx] = basis.nonredundant_indices[i]
        matrix.upper_to_mult[row_idx] = hashtable_insert!(ht, etmp)
        symbol_ht.labels[uprows[row_idx][1]] = UNKNOWN_PIVOT_COLUMN
    end
    trace.nonredundant_indices_before_reduce = basis.nonredundant_indices[1:(basis.n_nonredundant)]

    # needed for correct column count in symbol hashtable
    matrix.ncols_left = matrix.nrows_filled_upper

    f4_symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.labels[i] = UNKNOWN_PIVOT_COLUMN
    end

    matrix_fill_column_to_monom_map!(matrix, symbol_ht)

    linalg_autoreduce!(matrix, basis, params, trace, linalg=LinearAlgebra(:learn, :sparse))

    matrix_convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht, params)

    basis.n_filled = matrix.npivots + basis.n_processed
    basis.n_processed = matrix.npivots

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    k = 0
    i = 1
    @label Letsgo
    @inbounds while i <= basis.n_processed
        @inbounds for j in 1:k
            if hashtable_monom_is_divisible(
                basis.monoms[basis.n_filled - i + 1][1],
                basis.monoms[basis.nonredundant_indices[j]][1],
                ht
            )
                i += 1
                @goto Letsgo
            end
        end
        k += 1
        basis.nonredundant_indices[k] = basis.n_filled - i + 1
        basis.divmasks[k] = ht.divmasks[basis.monoms[basis.nonredundant_indices[k]][1]]
        i += 1
    end
    basis.n_nonredundant = k

    trace.output_nonredundant_indices = copy(basis.nonredundant_indices[1:k])
end

function f4_learn!(
    trace::Trace,
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @invariant basis_well_formed(ring, basis, hashtable)

    @invariant basis == trace.gb_basis
    @invariant params.reduced

    basis_make_monic!(basis, params.arithmetic, params.changematrix)

    matrix = matrix_initialize(ring, C)
    update_ht = hashtable_initialize_secondary(hashtable)
    symbol_ht = hashtable_initialize_secondary(hashtable)

    pairset_size = f4_update!(pairset, basis, hashtable, update_ht)

    while !isempty(pairset)
        degree_i, npairs_i = f4_select_critical_pairs!(pairset, basis, matrix, hashtable, symbol_ht)
        push!(trace.critical_pair_sequence, (degree_i, npairs_i))

        f4_symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)

        f4_reduction_learn!(trace, basis, matrix, hashtable, symbol_ht, params)

        pairset_size = f4_update!(pairset, basis, hashtable, update_ht)

        matrix_reinitialize!(matrix, 0)
        hashtable_reinitialize!(symbol_ht)
    end

    basis_move_redundant_elements!(basis)

    if params.reduced
        f4_reducegb_learn!(trace, ring, basis, matrix, hashtable, symbol_ht, params)
    end

    trace_finalize!(trace)

    perm = basis_standardize!(ring, basis, hashtable, params.arithmetic, false)
    trace.output_sort_indices = perm

    @invariant basis_well_formed(ring, basis, hashtable)

    nothing
end

###
# F4 apply stage

function matrix_fill_column_to_monom_map!(
    trace::Trace,
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable
)
    # monoms from symbolic table represent one column in the matrix
    load = symbol_ht.load

    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        if symbol_ht.labels[i] == PIVOT_COLUMN
            k += 1
        end
    end

    column_to_monom = matrix.column_to_monom

    matrix.ncols_left = k  # CHECK!
    # -1 as long as hashtable load is always 1 more than actual
    matrix.ncols_right = load - matrix.ncols_left - 1

    # store the other direction of mapping,
    # hash -> column
    @inbounds for k in 1:length(column_to_monom)
        symbol_ht.labels[column_to_monom[k]] = k
    end

    @inbounds for k in 1:(matrix.nrows_filled_upper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = symbol_ht.labels[row[j]]
        end
    end

    @inbounds for k in 1:(matrix.nrows_filled_lower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = symbol_ht.labels[row[j]]
        end
    end
end

function f4_reduction_apply!(
    trace::Trace,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    f4_iteration::Int,
    cache_column_order::Bool,
    params::AlgorithmParameters
)
    # We bypass the sorting of the matrix columns when apply is used the second
    # time
    if cache_column_order
        if length(trace.matrix_sorted_columns) >= f4_iteration
            matrix.column_to_monom = trace.matrix_sorted_columns[f4_iteration]
            matrix_fill_column_to_monom_map!(trace, matrix, symbol_ht)
        else
            matrix_fill_column_to_monom_map!(matrix, symbol_ht)
            push!(trace.matrix_sorted_columns, matrix.column_to_monom)
        end
    else
        matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    end

    flag = linalg_main!(matrix, basis, params, trace, linalg=LinearAlgebra(:apply, :sparse))
    !flag && return (false, false)

    matrix_convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht, params)

    # Check that the leading terms were not reduced to zero accidentally
    pivot_indices = trace.matrix_pivot_indices[f4_iteration]
    @inbounds for i in 1:(matrix.npivots)
        sgn = basis.monoms[basis.n_processed + i][1]
        sgn != pivot_indices[i] && return (false, false)
    end

    if cache_column_order
        matrix_pivot_signature =
            matrix_compute_pivot_signature(basis.monoms, basis.n_processed + 1, matrix.npivots)
        if matrix_pivot_signature != trace.matrix_pivot_signatures[f4_iteration]
            return false, false
        end
    end

    true, cache_column_order
end

function f4_symbolic_preprocessing!(
    trace::Trace,
    f4_iteration::Int,
    basis::Basis,
    matrix::MacaulayMatrix,
    hashtable::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    lowrows, lowmults = trace.matrix_lower_rows[f4_iteration]
    uprows, upmults = trace.matrix_upper_rows[f4_iteration]
    matrix_info = trace.matrix_infos[f4_iteration]
    nonzeroed_rows = trace.matrix_nonzeroed_rows[f4_iteration]
    nlow = length(nonzeroed_rows)
    nup = length(uprows)

    resize!(matrix.upper_rows, nup)
    resize!(matrix.lower_rows, nlow)
    resize!(matrix.lower_to_coeffs, nlow)
    resize!(matrix.upper_to_coeffs, nup)

    hashtable_resize_if_needed!(symbol_ht, nlow + nup + 2)
    @inbounds for i in 1:nlow
        mult_idx = lowmults[i]
        poly_idx = lowrows[i]

        h = hashtable.hashvals[mult_idx]
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        matrix.lower_rows[i] =
            matrix_polynomial_multiple_to_row!(matrix, symbol_ht, hashtable, h, etmp, rpoly)

        symbol_ht.labels[matrix.lower_rows[i][1]] = PIVOT_COLUMN

        matrix.lower_to_coeffs[i] = poly_idx
    end

    @invariant map(i -> length(matrix.lower_rows[i]), 1:nlow) ==
               map(i -> length(basis.coeffs[matrix.lower_to_coeffs[i]]), 1:nlow)

    @inbounds for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashvals[mult_idx]
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        matrix.upper_rows[i] =
            matrix_polynomial_multiple_to_row!(matrix, symbol_ht, hashtable, h, etmp, rpoly)

        symbol_ht.labels[matrix.upper_rows[i][1]] = PIVOT_COLUMN

        matrix.upper_to_coeffs[i] = poly_idx
    end

    @invariant map(i -> length(matrix.upper_rows[i]), 1:nup) ==
               map(i -> length(basis.coeffs[matrix.upper_to_coeffs[i]]), 1:nup)

    i = MonomId(symbol_ht.offset)
    @inbounds while i <= symbol_ht.load
        if symbol_ht.labels[i] == NON_PIVOT_COLUMN
            symbol_ht.labels[i] = UNKNOWN_PIVOT_COLUMN
        end
        i += MonomId(1)
    end

    matrix.nrows_filled_lower = nlow
    matrix.nrows_filled_upper = nup
end

function f4_autoreduce_apply!(
    trace::Trace,
    basis::Basis,
    matrix::MacaulayMatrix,
    hashtable::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    f4_iteration::Int,
    cache_column_order::Bool,
    params::AlgorithmParameters
) where {M}
    lowrows, _ = trace.matrix_lower_rows[end]
    uprows, upmults = trace.matrix_upper_rows[end]
    nlow = length(lowrows)
    nup = length(uprows)

    resize!(matrix.upper_rows, nup)
    resize!(matrix.lower_rows, nlow)
    resize!(matrix.lower_to_coeffs, nlow)
    resize!(matrix.upper_to_coeffs, nup)

    matrix.ncols_left = 0
    matrix.ncols_right = 0
    matrix.nrows_filled_upper = nup
    matrix.nrows_filled_lower = nlow

    etmp = monom_construct_const(M, hashtable.nvars)
    # etmp is now set to zero, and has zero hash

    hashtable_resize_if_needed!(symbol_ht, nlow + nup + 2)

    # needed for correct column count in symbol hashtable
    matrix.ncols_left = matrix.nrows_filled_upper

    @inbounds for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashvals[mult_idx]
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        # vidx = hashtable_insert!(symbol_ht, etmp)
        matrix.upper_rows[i] =
            matrix_polynomial_multiple_to_row!(matrix, symbol_ht, hashtable, h, etmp, rpoly)

        matrix.upper_to_coeffs[i] = poly_idx
    end

    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.labels[i] = UNKNOWN_PIVOT_COLUMN
    end

    matrix.nrows_filled_lower = nlow
    matrix.nrows_filled_upper = nup

    if cache_column_order
        if length(trace.matrix_sorted_columns) >= f4_iteration
            matrix.column_to_monom = trace.matrix_sorted_columns[f4_iteration]
            matrix_fill_column_to_monom_map!(trace, matrix, symbol_ht)
        else
            matrix_fill_column_to_monom_map!(matrix, symbol_ht)
            push!(trace.matrix_sorted_columns, matrix.column_to_monom)
        end
    else
        matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    end

    flag = linalg_autoreduce!(matrix, basis, params, linalg=LinearAlgebra(:apply, :sparse))
    !flag && return false

    matrix_convert_rows_to_basis_elements!(matrix, basis, hashtable, symbol_ht, params)

    basis.n_filled = matrix.npivots + basis.n_processed
    basis.n_processed = matrix.npivots

    output_nonredundant = trace.output_nonredundant_indices
    for i in 1:length(output_nonredundant)
        basis.nonredundant_indices[i] = output_nonredundant[i]
        basis.divmasks[i] = hashtable.divmasks[basis.monoms[basis.nonredundant_indices[i]][1]]
    end
    basis.n_nonredundant = length(output_nonredundant)
    true
end

function f4_standardize_basis_in_apply!(ring::PolyRing, trace::Trace, arithmetic)
    basis = trace.gb_basis
    buf = trace.buf_basis
    basis.n_processed = basis.n_filled = basis.n_nonredundant
    output_nonredundant = trace.output_nonredundant_indices
    output_sort = trace.output_sort_indices
    @inbounds for i in 1:(basis.n_processed)
        basis.coeffs[i] = buf.coeffs[output_nonredundant[output_sort[i]]]
    end
    buf.n_processed = buf.n_nonredundant = 0
    buf.n_filled = trace.input_basis.n_filled
    basis_make_monic!(basis, arithmetic, false)
end

function f4_apply!(
    trace::Trace,
    ring::PolyRing,
    basis::Basis{C},
    params::AlgorithmParameters
) where {C <: Coeff}
    ring = trace.ring
    @invariant basis_well_formed(ring, basis, trace.hashtable)
    @invariant params.reduced

    basis_make_monic!(basis, params.arithmetic, params.changematrix)

    iters_total = length(trace.matrix_infos) - 1
    iters = 0
    hashtable = trace.hashtable

    symbol_ht = hashtable_initialize_secondary(hashtable)
    matrix = matrix_initialize(ring, C)

    basis_update!(basis, hashtable)

    cache_column_order = true

    while iters < iters_total
        iters += 1

        f4_symbolic_preprocessing!(trace, iters, basis, matrix, hashtable, symbol_ht)

        flag, new_cache_column_order = f4_reduction_apply!(
            trace,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            iters,
            cache_column_order,
            params
        )
        if !flag
            # Unlucky cancellation of basis coefficients may have happened
            return false
        end
        cache_column_order = new_cache_column_order && cache_column_order

        basis_update!(basis, hashtable)

        hashtable_reinitialize!(symbol_ht)
    end

    if params.reduced
        symbol_ht = hashtable_initialize_secondary(hashtable)
        flag = f4_autoreduce_apply!(
            trace,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            iters_total + 1,
            cache_column_order,
            params
        )
        if !flag
            return false
        end
    end

    f4_standardize_basis_in_apply!(ring, trace, params.arithmetic)
    basis = trace.gb_basis

    @invariant basis_well_formed(ring, basis, hashtable)

    trace.napply += 1

    true
end
