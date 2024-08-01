# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# F4 learn & apply with tracing
#
# Tracing is implemented as an add-on to the F4 algorithm from f4.jl. Here, we
# try to reuse the functions from f4.jl as much as possible. At the same time,
# f4.jl is completely independent from this file.

###
# F4 learn stage

# Same as f4_initialize_structs, but also initializes a trace
function f4_initialize_structs_with_trace(
    ctx::Context,
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    make_monic=true,
    sort_input=true
) where {M <: Monom, C <: Coeff}
    basis, pairset, hashtable, permutation = f4_initialize_structs(
        ctx,
        ring,
        monoms,
        coeffs,
        params,
        make_monic=make_monic,
        sort_input=sort_input
    )

    trace = trace_initialize(
        ctx,
        ring,
        basis_deepcopy(ctx, basis),
        basis,
        hashtable,
        permutation,
        params
    )

    trace, basis, pairset, hashtable, permutation
end

function standardize_basis_in_learn!(
    trace::TraceF4,
    ring::PolyRing,
    basis::Basis,
    ht::MonomialHashtable,
    ord::Any,
    arithmetic
)
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
    perm = sort_polys_by_lead_increasing!(basis, ht, false, ord=ord)
    trace.output_sort_indices = perm
    basis_make_monic!(basis, arithmetic, false)
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

function reduction_learn!(
    ctx::Context,
    trace::TraceF4,
    basis::Basis,
    matrix::MacaulayMatrix,
    hashtable::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters
)
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    linalg_main!(ctx, matrix, basis, params, trace, linalg=LinearAlgebra(:learn, :sparse))
    matrix_convert_rows_to_basis_elements!(
        ctx,
        matrix,
        basis,
        hashtable,
        symbol_ht,
        params;
        batched_ht_insert=true
    )
    pivot_indices =
        map(i -> Int32(basis.monoms[basis.nprocessed + i][1]), 1:(matrix.npivots))
    push!(trace.matrix_pivot_indices, pivot_indices)
    matrix_pivot_signature =
        matrix_compute_pivot_signature(basis.monoms, basis.nprocessed + 1, matrix.npivots)
    push!(trace.matrix_pivot_signatures, matrix_pivot_signature)
end

function f4_reducegb_learn!(
    ctx::Context,
    trace::TraceF4,
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom}
    etmp = monom_construct_const(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    matrix_reinitialize!(matrix, basis.nnonredundant)
    uprows = matrix.upper_rows

    # add all non redundant elements from basis
    # as matrix upper rows
    @inbounds for i in 1:(basis.nnonredundant) #
        row_idx = matrix.nrows_filled_upper += 1
        uprows[row_idx] = matrix_polynomial_multiple_to_row!(
            ctx,
            matrix,
            symbol_ht,
            ht,
            MonomHash(0),
            etmp,
            basis.monoms[basis.nonredundant[i]]
        )

        matrix.upper_to_coeffs[row_idx] = basis.nonredundant[i]
        matrix.upper_to_mult[row_idx] = hashtable_insert!(ht, etmp)
        # set lead index as 1
        hv = symbol_ht.hashdata[uprows[row_idx][1]]
        symbol_ht.hashdata[uprows[row_idx][1]] =
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask)
    end
    trace.nonredundant_indices_before_reduce = basis.nonredundant[1:(basis.nnonredundant)]

    # needed for correct column count in symbol hashtable
    matrix.ncols_left = matrix.nrows_filled_upper

    # @info "" matrix.nrows_filled_upper

    f4_symbolic_preprocessing!(ctx, basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        hv = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] = Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask)
    end

    matrix_fill_column_to_monom_map!(matrix, symbol_ht)

    # @info "" matrix.nrows_filled_upper

    linalg_autoreduce!(
        ctx,
        matrix,
        basis,
        params,
        trace,
        linalg=LinearAlgebra(:learn, :sparse)
    )

    matrix_convert_rows_to_basis_elements!(ctx, matrix, basis, ht, symbol_ht, params)

    # println(basis.nnonredundant)
    # println(basis.nprocessed)
    # println(basis.nfilled)

    # @info "" length(trace.matrix_upper_rows[end][1])

    basis.nfilled = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    k = 0
    i = 1
    @label Letsgo
    @inbounds while i <= basis.nprocessed
        @inbounds for j in 1:k
            if hashtable_monom_is_divisible(
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

    # println(basis.nnonredundant)

    trace.output_nonredundant_indices = copy(basis.nonredundant[1:k])
end

@timeit function f4_learn!(
    ctx::Context,
    trace::TraceF4,
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    # @invariant hashtable_well_formed(:input_f4!, ring, hashtable)
    @invariant basis_well_formed(ring, basis, hashtable)
    # @invariant pairset_well_formed(:input_f4!, pairset, basis, ht)

    @invariant basis == trace.gb_basis
    @invariant params.reduced

    @log :debug "Entering F4 Learn phase."
    basis_make_monic!(basis, params.arithmetic, params.changematrix)

    matrix = matrix_initialize(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = hashtable_initialize_secondary(hashtable)
    symbol_ht = hashtable_initialize_secondary(hashtable)

    # add the first batch of critical pairs to the pairset
    @log :debug "Processing initial polynomials, generating first critical pairs"
    pairset_size = f4_update!(ctx, pairset, basis, hashtable, update_ht)
    @log :debug "Out of $(basis.nfilled) polynomials, $(basis.nprocessed) are non-redundant"
    @log :debug "Generated $(pairset.load) critical pairs"

    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log :debug "F4: iteration $i"
        @log :debug "F4: available $(pairset.load) pairs"

        degree_i, npairs_i = f4_select_critical_pairs!(
            ctx,
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            0,
            maxpairs=params.maxpairs
        )
        push!(trace.critical_pair_sequence, (degree_i, npairs_i))

        f4_symbolic_preprocessing!(ctx, basis, matrix, hashtable, symbol_ht)

        # reduces polys and obtains new potential basis elements
        reduction_learn!(ctx, trace, basis, matrix, hashtable, symbol_ht, params)

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = f4_update!(ctx, pairset, basis, hashtable, update_ht)

        # clear symbolic hashtable
        # clear matrix
        # matrix    = matrix_initialize(ring, C)
        matrix_reinitialize!(matrix, 0)
        hashtable_reinitialize!(symbol_ht)
        # symbol_ht = hashtable_initialize_secondary(hashtable)

        if i > 10_000
            @log :warn "Something has gone wrong in the F4 learn stage. Error will follow."
            @log_memory_locals
            __throw_maximum_iterations_exceeded_in_f4(i)
        end
    end

    if params.sweep
        @log :debug "Sweeping redundant elements in the basis"
        basis_sweep_redundant!(basis, hashtable)
    end

    # mark redundant elements
    basis_mark_redundant_elements!(basis)
    @log :debug "Filtered elements marked redundant"

    if params.reduced
        @log :debug "Autoreducing the final basis.."
        f4_reducegb_learn!(ctx, trace, ring, basis, matrix, hashtable, symbol_ht, params)
    end

    @log :debug "Finalizing computation trace"
    trace_finalize!(trace)

    standardize_basis_in_learn!(
        trace,
        ring,
        basis,
        hashtable,
        hashtable.ord,
        params.arithmetic
    )

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    @invariant basis_well_formed(ring, basis, hashtable)

    nothing
end

###
# F4 apply stage

function matrix_fill_column_to_monom_map!(
    trace::TraceF4,
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable
)
    # monoms from symbolic table represent one column in the matrix
    hdata = symbol_ht.hashdata
    load = symbol_ht.load

    # number of pivotal cols
    k = 0
    @inbounds for i in (symbol_ht.offset):load
        if hdata[i].idx == PIVOT_COLUMN
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
        hv = hdata[column_to_monom[k]]
        hdata[column_to_monom[k]] = Hashvalue(k, hv.hash, hv.divmask)
    end

    @inbounds for k in 1:(matrix.nrows_filled_upper)
        row = matrix.upper_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end

    @inbounds for k in 1:(matrix.nrows_filled_lower)
        row = matrix.lower_rows[k]
        for j in 1:length(row)
            row[j] = hdata[row[j]].idx
        end
    end
end

function f4_reduction_apply!(
    ctx::Context,
    trace::TraceF4,
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
            # TODO: see "TODO: (I)" in src/groebner/groebner.jl
            matrix.column_to_monom = trace.matrix_sorted_columns[f4_iteration]
            matrix_fill_column_to_monom_map!(trace, matrix, symbol_ht)
        else
            matrix_fill_column_to_monom_map!(matrix, symbol_ht)
            push!(trace.matrix_sorted_columns, matrix.column_to_monom)
        end
    else
        matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    end

    flag = linalg_main!(
        ctx,
        matrix,
        basis,
        params,
        trace,
        linalg=LinearAlgebra(:apply, :sparse)
    )
    if !flag
        @log :warn "In apply, in linear algebra, some of the matrix rows reduced to zero."
        return (false, false)
    end

    matrix_convert_rows_to_basis_elements!(ctx, matrix, basis, ht, symbol_ht, params)

    # Check that the leading terms were not reduced to zero accidentally
    pivot_indices = trace.matrix_pivot_indices[f4_iteration]
    @inbounds for i in 1:(matrix.npivots)
        sgn = basis.monoms[basis.nprocessed + i][1]
        if sgn != pivot_indices[i]
            @log :debug "In apply, some leading terms cancelled out!"
            return (false, false)
        end
    end

    if cache_column_order
        matrix_pivot_signature = matrix_compute_pivot_signature(
            basis.monoms,
            basis.nprocessed + 1,
            matrix.npivots
        )
        if matrix_pivot_signature != trace.matrix_pivot_signatures[f4_iteration]
            @log :info """
            In apply, on iteration $(f4_iteration) of F4, some trailing terms cancelled out.
            hash (expected):    $(trace.matrix_pivot_signatures[f4_iteration])
            hash (got):         $(matrix_pivot_signature)"""
            # return false, false
            return false, false
        end
    end

    true, cache_column_order
end

function f4_symbolic_preprocessing!(
    ctx::Context,
    trace::TraceF4,
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
    nlow = length(nonzeroed_rows) # matrix_info.nlow # length(nonzeroed_rows)
    nup = length(uprows) # matrix_info.nup # length(uprows)

    resize!(matrix.upper_rows, nup)
    resize!(matrix.lower_rows, nlow)
    resize!(matrix.lower_to_coeffs, nlow)
    resize!(matrix.upper_to_coeffs, nup)
    # matrix.upper_rows = Vector{Vector{ColumnLabel}}(undef, nup)
    # matrix.lower_rows = Vector{Vector{ColumnLabel}}(undef, nlow)
    # matrix.lower_to_coeffs = Vector{Int}(undef, nlow)
    # matrix.upper_to_coeffs = Vector{Int}(undef, nup)

    hashtable_resize_if_needed!(symbol_ht, nlow + nup + 2)
    @inbounds for i in 1:nlow
        mult_idx = lowmults[i]
        poly_idx = lowrows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        matrix.lower_rows[i] = matrix_polynomial_multiple_to_row!(
            ctx,
            matrix,
            symbol_ht,
            hashtable,
            h,
            etmp,
            rpoly
        )

        hv = symbol_ht.hashdata[matrix.lower_rows[i][1]]
        symbol_ht.hashdata[matrix.lower_rows[i][1]] =
            Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask)

        matrix.lower_to_coeffs[i] = poly_idx
    end

    @invariant map(i -> length(matrix.lower_rows[i]), 1:nlow) ==
               map(i -> length(basis.coeffs[matrix.lower_to_coeffs[i]]), 1:nlow)

    @inbounds for i in 1:nup
        mult_idx = upmults[i]
        poly_idx = uprows[i]

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        matrix.upper_rows[i] = matrix_polynomial_multiple_to_row!(
            ctx,
            matrix,
            symbol_ht,
            hashtable,
            h,
            etmp,
            rpoly
        )

        hv = symbol_ht.hashdata[matrix.upper_rows[i][1]]
        symbol_ht.hashdata[matrix.upper_rows[i][1]] =
            Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask)

        matrix.upper_to_coeffs[i] = poly_idx
    end

    @invariant map(i -> length(matrix.upper_rows[i]), 1:nup) ==
               map(i -> length(basis.coeffs[matrix.upper_to_coeffs[i]]), 1:nup)

    i = MonomId(symbol_ht.offset)
    @inbounds while i <= symbol_ht.load
        if symbol_ht.hashdata[i].idx == NON_PIVOT_COLUMN
            hv = symbol_ht.hashdata[i]
            symbol_ht.hashdata[i] = Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask)
        end
        i += MonomId(1)
    end

    matrix.nrows_filled_lower = nlow
    matrix.nrows_filled_upper = nup
end

function f4_autoreduce_apply!(
    ctx::Context,
    trace::TraceF4,
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
    # nlow = trace.matrix_infos[end].nlow
    # nup = trace.matrix_infos[end].nup

    # matrix.upper_rows = Vector{Vector{ColumnLabel}}(undef, nup)
    # matrix.lower_rows = Vector{Vector{ColumnLabel}}(undef, nlow)
    # matrix.lower_to_coeffs = Vector{Int}(undef, nlow)
    # matrix.upper_to_coeffs = Vector{Int}(undef, nup)
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

        h = hashtable.hashdata[mult_idx].hash
        etmp = hashtable.monoms[mult_idx]
        rpoly = basis.monoms[poly_idx]

        # vidx = hashtable_insert!(symbol_ht, etmp)
        matrix.upper_rows[i] = matrix_polynomial_multiple_to_row!(
            ctx,
            matrix,
            symbol_ht,
            hashtable,
            h,
            etmp,
            rpoly
        )

        matrix.upper_to_coeffs[i] = poly_idx
    end

    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        hv = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] = Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask)
    end

    matrix.nrows_filled_lower = nlow
    matrix.nrows_filled_upper = nup

    if cache_column_order
        if length(trace.matrix_sorted_columns) >= f4_iteration
            # TODO: see "TODO: (I)" in src/groebner/groebner.jl
            matrix.column_to_monom = trace.matrix_sorted_columns[f4_iteration]
            matrix_fill_column_to_monom_map!(trace, matrix, symbol_ht)
        else
            matrix_fill_column_to_monom_map!(matrix, symbol_ht)
            push!(trace.matrix_sorted_columns, matrix.column_to_monom)
        end
    else
        matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    end

    flag = linalg_autoreduce!(
        ctx,
        matrix,
        basis,
        params,
        linalg=LinearAlgebra(:apply, :sparse)
    )
    if !flag
        @log :warn "In apply, the final autoreduction of the basis failed"
        return false
    end
    matrix_convert_rows_to_basis_elements!(ctx, matrix, basis, hashtable, symbol_ht, params)

    basis.nfilled = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    output_nonredundant = trace.output_nonredundant_indices
    for i in 1:length(output_nonredundant)
        basis.nonredundant[i] = output_nonredundant[i]
        basis.divmasks[i] =
            hashtable.hashdata[basis.monoms[basis.nonredundant[i]][1]].divmask
    end
    basis.nnonredundant = length(output_nonredundant)
    true
end

function f4_standardize_basis_in_apply!(ring::PolyRing, trace::TraceF4, arithmetic)
    basis = trace.gb_basis
    buf = trace.buf_basis
    basis.size = basis.nprocessed = basis.nfilled = basis.nnonredundant
    output_nonredundant = trace.output_nonredundant_indices
    output_sort = trace.output_sort_indices
    @inbounds for i in 1:(basis.nprocessed)
        basis.coeffs[i] = buf.coeffs[output_nonredundant[output_sort[i]]]
    end
    buf.nprocessed = buf.nnonredundant = 0
    buf.nfilled = trace.input_basis.nfilled
    basis_make_monic!(basis, arithmetic, false)
end

@timeit function f4_apply!(
    ctx::Context,
    trace::TraceF4,
    ring::PolyRing,
    basis::Basis{C},
    params::AlgorithmParameters
) where {C <: Coeff}
    @invariant basis_well_formed(ring, basis, trace.hashtable)
    @invariant params.reduced

    basis_make_monic!(basis, params.arithmetic, params.changematrix)

    iters_total = length(trace.matrix_infos) - 1
    iters = 0
    hashtable = trace.hashtable

    symbol_ht = hashtable_initialize_secondary(hashtable)
    matrix = matrix_initialize(ring, C)
    @log :debug "Applying F4 modulo $(ring.ch)"
    @log :debug "Using parameters" params.arithmetic
    @invariant (
        _T = typeof(divisor(params.arithmetic)); _T(ring.ch) == divisor(params.arithmetic)
    )

    basis_update!(basis, hashtable)

    cache_column_order = true

    while iters < iters_total
        iters += 1
        @log :debug "F4 Apply iteration $iters"

        f4_symbolic_preprocessing!(ctx, trace, iters, basis, matrix, hashtable, symbol_ht)

        flag, new_cache_column_order = f4_reduction_apply!(
            ctx,
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

        @log_memory_locals
    end

    # basis_mark_redundant_elements!(basis)

    if params.reduced
        @log :debug "Autoreducing the final basis.."
        symbol_ht = hashtable_initialize_secondary(hashtable)
        flag = f4_autoreduce_apply!(
            ctx,
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

    true
end
