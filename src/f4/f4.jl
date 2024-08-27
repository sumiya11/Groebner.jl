# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

### 
# Main file that defines the f4! function.

function f4_initialize_structs(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    make_monic=true,
    sort_input=true
) where {M <: Monom, C <: Coeff}
    tablesize = hashtable_select_initial_size(ring, monoms)
    basis = basis_initialize(ring, length(monoms), C)
    pairset = pairset_initialize(monom_entrytype(M))
    hashtable = hashtable_initialize(ring, params.rng, M, tablesize)

    basis_fill_data!(basis, hashtable, monoms, coeffs)
    hashtable_fill_divmasks!(hashtable)

    if params.changematrix
        basis_changematrix_initialize!(basis, hashtable)
    end

    if sort_input
        # The sorting of input polynomials is not deterministic across different
        # Julia versions when sorting only w.r.t. the leading term.
        permutation = sort_polys_by_lead_increasing!(basis, hashtable, params.changematrix)
    else
        permutation = collect(1:(basis.nfilled))
    end

    # We do not need monic polynomials when computing normal forms
    if make_monic
        basis_make_monic!(basis, params.arithmetic, params.changematrix)
    end

    basis, pairset, hashtable, permutation
end

function f4_reduction!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters
)
    if SEMIGROUP_ON[]
        d = symbol_ht.monoms[1][1]
        for i in 2:(symbol_ht.load)
            @invariant semigroup_check_normalized_monom(symbol_ht.monoms[i])
            @invariant symbol_ht.monoms[i][1] == d
        end
    end

    matrix_fill_column_to_monom_map!(matrix, symbol_ht)

    push!(DATA[:matrix_size], size(matrix))

    linalg_main!(matrix, basis, params)

    push!(DATA[:useful_rows], round(matrix.npivots / matrix.nrows_filled_lower, digits=10))

    matrix_convert_rows_to_basis_elements!(
        matrix,
        basis,
        ht,
        symbol_ht,
        params;
        batched_ht_insert=true
    )
end

function f4_update!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable,
    update_ht::MonomialHashtable
)
    @invariant basis.nfilled >= basis.nprocessed
    @inbounds for i in (basis.nprocessed + 1):(basis.nfilled)
        basis_is_new_polynomial_redundant!(pairset, basis, ht, update_ht, i) && continue
        pairset_resize_lcms_if_needed!(pairset, basis.nfilled)
        pairset_resize_if_needed!(pairset, basis.nfilled)
        pairset_update!(pairset, basis, ht, update_ht, i)
    end
    basis_update!(basis, ht)
end

function f4_symbolic_preprocessing!(
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    arithmetic
)
    # Monomials that represent the columns of the matrix are stored in the
    # symbol_ht hashtable.
    symbol_load = symbol_ht.load
    ncols = matrix.ncols_left
    matrix_resize_upper_part_if_needed!(matrix, ncols + symbol_load)

    if SEMIGROUP_ON[]
        for i in 2:(symbol_ht.load)
            @invariant semigroup_check_normalized_monom(symbol_ht.monoms[i])
        end
    end

    # Traverse all monomials in symbol_ht and search for a polynomial reducer
    # for each monomial. The hashtable grows as polynomials with new monomials
    # are added to the matrix, and the loop accounts for that.
    i = symbol_ht.offset
    @inbounds while i <= symbol_ht.load
        if symbol_ht.hashdata[i].idx != NON_PIVOT_COLUMN
            i += 1
            continue
        end
        matrix_resize_upper_part_if_needed!(matrix, matrix.nrows_filled_upper + 1)

        hashval = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] =
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hashval.hash, hashval.divmask)
        matrix.ncols_left += 1
        f4_find_multiplied_reducer!(basis, matrix, ht, symbol_ht, MonomId(i), arithmetic)
        i += 1
    end

    nothing
end

function f4_autoreduce!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    params
) where {M <: Monom}
    etmp = monom_construct_const(M, ht.nvars)

    matrix_reinitialize!(matrix, basis.nnonredundant)
    uprows = matrix.upper_rows
    lowrows = matrix.lower_rows

    @inbounds for i in 1:(basis.nnonredundant)
        if SEMIGROUP_ON[]
            if semigroup_is_a_relation_lead(
                ht.monoms[basis.monoms[basis.nonredundant[i]][1]]
            )
                continue
            end
        end
        matrix.nrows_filled_upper += 1
        row_idx = matrix.nrows_filled_upper
        uprows[row_idx] = matrix_polynomial_multiple_to_row!(
            matrix,
            symbol_ht,
            ht,
            MonomHash(0),
            etmp,
            basis.monoms[basis.nonredundant[i]]
        )
        matrix.upper_to_coeffs[row_idx] = basis.nonredundant[i]
        matrix.upper_to_mult[row_idx] = hashtable_insert!(ht, etmp)
        hv = symbol_ht.hashdata[uprows[row_idx][1]]
        symbol_ht.hashdata[uprows[row_idx][1]] =
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask)
    end

    f4_symbolic_preprocessing!(basis, matrix, ht, symbol_ht, params.arithmetic)

    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        hv = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] = Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask)
    end

    matrix_fill_column_to_monom_map!(matrix, symbol_ht)

    if SEMIGROUP_ON[]
        for i in 2:(symbol_ht.load)
            @invariant semigroup_check_normalized_monom(symbol_ht.monoms[i])
        end
    end

    linalg_autoreduce!(matrix, basis, params)

    matrix_convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht, params)

    basis.nfilled = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

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
end

function f4_select_tobereduced!(
    basis::Basis,
    tobereduced::Basis{C},
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable{M},
    ht::MonomialHashtable{M}
) where {C, M}
    matrix_reinitialize!(matrix, max(basis.nfilled, tobereduced.nfilled))
    resize!(matrix.lower_rows, tobereduced.nfilled)
    resize!(matrix.some_coeffs, tobereduced.nfilled)

    etmp = monom_construct_const(M, ht.nvars)

    @inbounds for i in 1:(tobereduced.nfilled)
        matrix.nrows_filled_lower += 1
        row_idx = matrix.nrows_filled_lower

        gen = tobereduced.monoms[i]
        h = MonomHash(0)
        matrix.lower_rows[row_idx] =
            matrix_polynomial_multiple_to_row!(matrix, symbol_ht, ht, h, etmp, gen)
        matrix.lower_to_coeffs[row_idx] = i
        # TODO: not really needed here
        matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)
        matrix.some_coeffs[row_idx] = tobereduced.coeffs[i]
    end

    basis.nnonredundant = basis.nprocessed = basis.nfilled
    basis.isredundant .= 0
    @inbounds for i in 1:(basis.nnonredundant)
        basis.nonredundant[i] = i
        basis.divmasks[i] = ht.hashdata[basis.monoms[i][1]].divmask
    end

    nothing
end

function f4_find_lead_divisor_use_divmask(i, divmask, basis)
    lead_divmasks = basis.divmasks
    @inbounds while i <= basis.nnonredundant
        if divmask_is_probably_divisible(divmask, lead_divmasks[i])
            break
        end
        i += 1
    end
    i
end

function f4_find_lead_divisor(i, monom, basis, ht)
    @inbounds while i <= basis.nnonredundant
        lead_monom = ht.monoms[basis.monoms[basis.nonredundant[i]][1]]
        if monom_is_divisible(monom, lead_monom)
            break
        end
        i += 1
    end
    i
end

# Finds a polynomial from the basis with the leading term that divides the given
# monomial. If such a polynomial is found, writes a multiple of it to the
# symbolic hashtable.
function f4_find_multiplied_reducer!(
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    monomid::MonomId,
    arithmetic
)
    @inbounds monom = symbol_ht.monoms[monomid]
    @inbounds quotient = ht.monoms[1]
    @inbounds divmask = symbol_ht.hashdata[monomid].divmask

    # Start of the search loop.
    i = 1
    @label Letsgo

    if ht.use_divmask
        i = f4_find_lead_divisor_use_divmask(i, divmask, basis)
    else
        i = f4_find_lead_divisor(i, monom, basis, ht)
    end

    # Reducer is not found, yield.
    i > basis.nnonredundant && return nothing

    # Here, we have found a polynomial from the basis with the leading monom
    # that divides the given monom.
    @inbounds poly = basis.monoms[basis.nonredundant[i]]
    @inbounds lead = ht.monoms[poly[1]]

    success, quotient = monom_is_divisible!(quotient, monom, lead)
    if !success # division mask failed for some reason
        i += 1
        @goto Letsgo
    end

    if SEMIGROUP_ON[]
        if semigroup_is_a_relation_lead(lead)
            i += 1
            @goto Letsgo
        end
    end

    # Using the fact that the hash is linear.
    @inbounds quotient_hash = symbol_ht.hashdata[monomid].hash - ht.hashdata[poly[1]].hash
    @invariant quotient_hash == monom_hash(quotient, ht.hasher)

    hashtable_resize_if_needed!(ht, length(poly))
    row = matrix_polynomial_multiple_to_row!(
        matrix,
        symbol_ht,
        ht,
        quotient_hash,
        quotient,
        poly
    )
    matrix.nrows_filled_upper += 1
    row_id = matrix.nrows_filled_upper
    @inbounds matrix.upper_rows[row_id] = row
    @inbounds matrix.upper_to_coeffs[row_id] = basis.nonredundant[i]

    if SEMIGROUP_ON[]
        sup = matrix.upper_rows[row_id]
        cfs = basis.coeffs[basis.nonredundant[i]]
        if !issorted(
            matrix.upper_rows[row_id],
            lt=(i, j) ->
                monom_isless(symbol_ht.monoms[i], symbol_ht.monoms[j], symbol_ht.ord),
            rev=true
        )
            perm = collect(1:length(matrix.upper_rows[row_id]))
            # sort!(
            #     perm,
            #     lt=(i, j) -> monom_isless(
            #         symbol_ht.monoms[matrix.upper_rows[row_id][i]],
            #         symbol_ht.monoms[matrix.upper_rows[row_id][j]],
            #         symbol_ht.ord
            #     ),
            #     rev=true
            # )
            sup = matrix.upper_rows[row_id][perm]
            cfs = basis.coeffs[basis.nonredundant[i]][perm]
            linalg_row_make_monic!(cfs, arithmetic)
        end
        @assert sup[1] == monomid
        @assert cfs[1] == 1
        matrix.upper_rows[row_id] = sup
        matrix.semigroup_upper_coeffs[row_id] = cfs
        matrix.upper_to_coeffs[row_id] = row_id
    end

    # Insert the quotient into the main hashtable.
    # TODO: This line is here with one sole purpose -- to support tracing.
    #       Probably want to factor it out.
    @inbounds matrix.upper_to_mult[row_id] = hashtable_insert!(ht, quotient)

    # Mark the current monomial as a pivot since a reducer has been found.
    @inbounds monom_hashval = symbol_ht.hashdata[monomid]
    @inbounds symbol_ht.hashdata[monomid] =
        Hashvalue(PIVOT_COLUMN, monom_hashval.hash, monom_hashval.divmask)

    nothing
end

function f4_select_critical_pairs!(
    pairset::Pairset{ExponentType},
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    arithmetic;
    maxpairs::Int=INT_INF,
    select_all::Bool=false
) where {ExponentType}
    # TODO: Why is this code type unstable !???
    npairs::Int = pairset.load
    if !select_all
        npairs = pairset_lowest_degree_pairs!(pairset)
    end
    npairs = min(npairs, maxpairs)

    if SEMIGROUP_ON[]
        # if npairs > 250
        #     @info "" npairs
        #     npairs = max(1, ceil(Int, npairs * 0.2))
        #     @info "" npairs
        # end
    end

    @invariant npairs > 0

    ps = pairset.pairs
    degs = pairset.degrees
    @inbounds deg::ExponentType = degs[1]

    sort_pairset_by_lcm!(pairset, npairs, ht)

    # When there is a limit on the number of selected pairs, we still add pairs
    # which have the same lcm as the selected ones.
    if npairs > maxpairs
        navailable = npairs
        npairs = maxpairs
        lastlcm = ps[npairs].lcm
        while npairs < navailable && ps[npairs + 1].lcm == lastlcm
            npairs += 1
        end
    end

    push!(DATA[:degree], deg)

    f4_add_critical_pairs_to_matrix!(
        pairset,
        npairs,
        basis,
        matrix,
        ht,
        symbol_ht,
        arithmetic
    )

    if SEMIGROUP_ON[]
        d = symbol_ht.monoms[1][1]
        for i in 2:(symbol_ht.load)
            @invariant semigroup_check_normalized_monom(symbol_ht.monoms[i])
            @invariant symbol_ht.monoms[i][1] == d
        end
    end

    # Remove selected parirs from the pairset.
    @inbounds for i in 1:(pairset.load - npairs)
        ps[i] = ps[i + npairs]
        degs[i] = degs[i + npairs]
    end
    pairset.load -= npairs

    push!(DATA[:pairs], npairs)
    
    deg, npairs
end

function f4_add_critical_pairs_to_matrix!(
    pairset::Pairset,
    npairs::Int,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    arithmetic
)
    matrix_reinitialize!(matrix, npairs)
    pairs = pairset.pairs
    uprows = matrix.upper_rows
    lowrows = matrix.lower_rows

    polys = Vector{Int}(undef, 2 * npairs)
    @inbounds etmp = ht.monoms[1]

    if SEMIGROUP_ON[]
        for i in 1:npairs
            lcm = pairs[i].lcm
            # lcm = hashtable_insert!(ht, semigroup_normalize_monom(ht.monoms[lcm]))
            p1, p2 = pairs[i].poly1, pairs[i].poly2

            if semigroup_is_a_relation_lead(ht.monoms[basis.monoms[p2][1]])
                p1, p2 = p2, p1
            end
            if semigroup_is_a_relation_lead(ht.monoms[basis.monoms[p2][1]])
                continue
            end

            if semigroup_is_a_relation_lead(ht.monoms[basis.monoms[p1][1]])
                vidx = basis.monoms[p2][1]
                etmp = monom_division!(etmp, ht.monoms[lcm], ht.monoms[vidx])
                htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

                matrix.nrows_filled_lower += 1
                row_idx = matrix.nrows_filled_lower
                lowrows[row_idx] = matrix_polynomial_multiple_to_row!(
                    matrix,
                    symbol_ht,
                    ht,
                    htmp,
                    etmp,
                    basis.monoms[p2]
                )
                matrix.lower_to_coeffs[row_idx] = p2
                matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)

                sup = lowrows[row_idx]
                cfs = basis.coeffs[p2]
                if !issorted(
                    lowrows[row_idx],
                    lt=(i, j) -> monom_isless(
                        symbol_ht.monoms[i],
                        symbol_ht.monoms[j],
                        symbol_ht.ord
                    ),
                    rev=true
                )
                    perm = collect(1:length(lowrows[row_idx]))
                    row = lowrows[row_idx]
                    # sort!(
                    #     perm,
                    #     lt=(i, j) -> @inbounds(monom_isless(
                    #         symbol_ht.monoms[row[i]],
                    #         symbol_ht.monoms[row[j]],
                    #         symbol_ht.ord
                    #     )),
                    #     rev=true
                    # )
                    sup = lowrows[row_idx][perm]
                    cfs = basis.coeffs[p2][perm]
                    linalg_row_make_monic!(cfs, arithmetic)
                end
                @assert cfs[1] == 1
                lowrows[row_idx] = sup
                matrix.semigroup_lower_coeffs[row_idx] = cfs
                matrix.lower_to_coeffs[row_idx] = row_idx

                @assert symbol_ht.hashdata[lowrows[row_idx][1]].idx == NON_PIVOT_COLUMN
            else
                vidx = basis.monoms[p1][1]
                etmp = monom_division!(etmp, ht.monoms[lcm], ht.monoms[vidx])
                htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

                matrix.nrows_filled_lower += 1
                row_idx = matrix.nrows_filled_lower
                lowrows[row_idx] = matrix_polynomial_multiple_to_row!(
                    matrix,
                    symbol_ht,
                    ht,
                    htmp,
                    etmp,
                    basis.monoms[p1]
                )
                matrix.lower_to_coeffs[row_idx] = p1
                matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)

                sup = lowrows[row_idx]
                cfs = basis.coeffs[p1]
                if !issorted(
                    lowrows[row_idx],
                    lt=(i, j) -> monom_isless(
                        symbol_ht.monoms[i],
                        symbol_ht.monoms[j],
                        symbol_ht.ord
                    ),
                    rev=true
                )
                    perm = collect(1:length(lowrows[row_idx]))
                    row = lowrows[row_idx]
                    # sort!(
                    #     perm,
                    #     lt=(i, j) -> @inbounds(monom_isless(
                    #         symbol_ht.monoms[row[i]],
                    #         symbol_ht.monoms[row[j]],
                    #         symbol_ht.ord
                    #     )),
                    #     rev=true
                    # )
                    sup = lowrows[row_idx][perm]
                    cfs = basis.coeffs[p1][perm]
                    linalg_row_make_monic!(cfs, arithmetic)
                end
                @assert cfs[1] == 1
                @assert allunique(sup)
                lowrows[row_idx] = sup
                matrix.semigroup_lower_coeffs[row_idx] = cfs
                matrix.lower_to_coeffs[row_idx] = row_idx

                vidx = basis.monoms[p2][1]
                etmp = monom_division!(etmp, ht.monoms[lcm], ht.monoms[vidx])
                htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash
                matrix.nrows_filled_lower += 1
                row_idx = matrix.nrows_filled_lower
                lowrows[row_idx] = matrix_polynomial_multiple_to_row!(
                    matrix,
                    symbol_ht,
                    ht,
                    htmp,
                    etmp,
                    basis.monoms[p2]
                )
                matrix.lower_to_coeffs[row_idx] = p2
                matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)

                sup = lowrows[row_idx]
                cfs = basis.coeffs[p2]
                if !issorted(
                    lowrows[row_idx],
                    lt=(i, j) -> monom_isless(
                        symbol_ht.monoms[i],
                        symbol_ht.monoms[j],
                        symbol_ht.ord
                    ),
                    rev=true
                )
                    perm = collect(1:length(lowrows[row_idx]))
                    # sort!(
                    #     perm,
                    #     lt=(i, j) -> monom_isless(
                    #         symbol_ht.monoms[lowrows[row_idx][i]],
                    #         symbol_ht.monoms[lowrows[row_idx][j]],
                    #         symbol_ht.ord
                    #     ),
                    #     rev=true
                    # )
                    sup = lowrows[row_idx][perm]
                    cfs = basis.coeffs[p2][perm]
                    linalg_row_make_monic!(cfs, arithmetic)
                end
                @assert cfs[1] == 1
                # @assert allunique(sup)
                lowrows[row_idx] = sup
                matrix.semigroup_lower_coeffs[row_idx] = cfs
                matrix.lower_to_coeffs[row_idx] = row_idx

                @assert symbol_ht.hashdata[lowrows[row_idx][1]].idx == NON_PIVOT_COLUMN
            end
        end

        return
    else
        for i in 1:npairs
            lcm = pairs[i].lcm
            p1, p2 = pairs[i].poly1, pairs[i].poly2
            
            vidx = basis.monoms[p1][1]
            etmp = monom_division!(etmp, ht.monoms[lcm], ht.monoms[vidx])
            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            matrix.nrows_filled_lower += 1
            row_idx = matrix.nrows_filled_lower
            lowrows[row_idx] = matrix_polynomial_multiple_to_row!(
                matrix,
                symbol_ht,
                ht,
                htmp,
                etmp,
                basis.monoms[p1]
            )
            matrix.lower_to_coeffs[row_idx] = p1
            matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)

            vidx = basis.monoms[p2][1]
            etmp = monom_division!(etmp, ht.monoms[lcm], ht.monoms[vidx])
            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash
            matrix.nrows_filled_lower += 1
            row_idx = matrix.nrows_filled_lower
            lowrows[row_idx] = matrix_polynomial_multiple_to_row!(
                matrix,
                symbol_ht,
                ht,
                htmp,
                etmp,
                basis.monoms[p2]
            )
            matrix.lower_to_coeffs[row_idx] = p2
            matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)
        end

        return
    end

    i = 1
    @inbounds while i <= npairs
        npolys = 1
        lcm = pairs[i].lcm
        j = i
        while j <= npairs && pairs[j].lcm == lcm
            polys[npolys] = pairs[j].poly1
            npolys += 1
            polys[npolys] = pairs[j].poly2
            npolys += 1
            j += 1
        end
        npolys -= 1

        sort_generators_by_position!(polys, npolys)

        prev = polys[1]
        vidx = basis.monoms[prev][1]
        etmp = monom_division!(etmp, ht.monoms[lcm], ht.monoms[vidx])
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        matrix.ncols_left += 1

        matrix.nrows_filled_upper += 1
        row_idx = matrix.nrows_filled_upper
        uprows[row_idx] = matrix_polynomial_multiple_to_row!(
            matrix,
            symbol_ht,
            ht,
            htmp,
            etmp,
            basis.monoms[prev]
        )
        matrix.upper_to_coeffs[row_idx] = prev
        matrix.upper_to_mult[row_idx] = hashtable_insert!(ht, etmp)

        hv = symbol_ht.hashdata[uprows[row_idx][1]]
        symbol_ht.hashdata[uprows[row_idx][1]] =
            Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask)

        for k in 1:npolys
            polys[k] == prev && continue

            # hashtable could have been reallocated, so refresh the pointer.
            elcm = ht.monoms[lcm]

            prev = polys[k]
            vidx = basis.monoms[prev][1]
            eidx = ht.monoms[vidx]
            etmp = monom_division!(etmp, elcm, eidx)
            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            matrix.nrows_filled_lower += 1
            row_idx = matrix.nrows_filled_lower
            lowrows[row_idx] = matrix_polynomial_multiple_to_row!(
                matrix,
                symbol_ht,
                ht,
                htmp,
                etmp,
                basis.monoms[prev]
            )
            matrix.lower_to_coeffs[row_idx] = prev
            matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)

            hv = symbol_ht.hashdata[lowrows[row_idx][1]]
            symbol_ht.hashdata[lowrows[row_idx][1]] =
                Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask)
        end

        i = j
    end
end

function f4!(
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @invariant basis_well_formed(ring, basis, hashtable)

    [empty!(v) for (k, v) in DATA]
    [TIME[k] = zero(v) for (k, v) in TIME]

    TIME[:total] = time_ns()

    basis_make_monic!(basis, params.arithmetic, params.changematrix)

    SEMIGROUP_ON[] && (@invariant semigroup_normalized(basis, hashtable))

    matrix = matrix_initialize(ring, C)
    update_ht = hashtable_initialize_secondary(hashtable)
    symbol_ht = hashtable_initialize_secondary(hashtable)

    f4_update!(pairset, basis, hashtable, update_ht)

    i = 1
    while !isempty(pairset)
        push!(DATA[:i], i)
        i += 1

        f4_select_critical_pairs!(
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            params.arithmetic,
            maxpairs=params.maxpairs
        )

        f4_symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht, params.arithmetic)
        
        TIME[:f4_reduction!] += @elapsed f4_reduction!(ring, basis, matrix, hashtable, symbol_ht, params)

        f4_update!(pairset, basis, hashtable, update_ht)

        matrix_reinitialize!(matrix, 0)
        hashtable_reinitialize!(symbol_ht)

        SEMIGROUP_ON[] && (@invariant semigroup_normalized(basis, hashtable))
    end

    if params.sweep
        basis_sweep_redundant!(basis, hashtable)
    end

    basis_mark_redundant_elements!(basis)

    SEMIGROUP_ON[] && (@invariant semigroup_normalized(basis, hashtable))

    if params.reduced
        f4_autoreduce!(ring, basis, matrix, hashtable, symbol_ht, params)
    end

    basis_standardize!(ring, basis, hashtable, params.arithmetic, params.changematrix)

    SEMIGROUP_ON[] && (@invariant semigroup_normalized(basis, hashtable))

    @invariant basis_well_formed(ring, basis, hashtable)

    TIME[:total] = ((time_ns() - TIME[:total])) / 1e9

    nothing
end

function f4_isgroebner!(
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable,
    arithmetic::AbstractArithmetic
) where {C <: Coeff}
    @invariant basis_well_formed(ring, basis, hashtable)
    basis_make_monic!(basis, arithmetic, false)
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(hashtable)
    update_ht = hashtable_initialize_secondary(hashtable)
    f4_update!(pairset, basis, hashtable, update_ht)
    isempty(pairset) && return true
    f4_select_critical_pairs!(
        pairset,
        basis,
        matrix,
        hashtable,
        symbol_ht,
        arithmetic,
        select_all=true
    )
    f4_symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht, arithmetic)
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    linalg_isgroebner!(matrix, basis, arithmetic)
end

function f4_normalform!(
    ring::PolyRing,
    basis::Basis{C},
    tobereduced::Basis{C},
    ht::MonomialHashtable,
    arithmetic::AbstractArithmetic
) where {C <: Coeff}
    @invariant basis_well_formed(ring, basis, ht)
    basis_make_monic!(basis, arithmetic, false)
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(ht)
    f4_select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)
    f4_symbolic_preprocessing!(basis, matrix, ht, symbol_ht, arithmetic)
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    linalg_normalform!(matrix, basis, arithmetic)
    matrix_convert_rows_to_basis_elements_nf!(matrix, tobereduced, ht, symbol_ht)
    nothing
end
