# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

### 
# Main file that defines the f4! function.

# Functions here mostly work with a subset of these objects:
# ring      - current polynomial ring,
# basis     - a struct that stores polynomials,
# matrix    - a struct that is used for F4-style reduction,
# hashtable - a hashtable that stores monomials.

@noinline __f4_throw_maximum_iterations_exceeded(iters) =
    throw("""Something probably went wrong in Groebner.jl/F4. 
          The number of F4 iterations exceeded $iters. 
          Please consider submitting a GitHub issue.""")

# Given the polynomial ring and the arrays of monomials and coefficients,
# initializes and returns the structs that are necessary for calling F4
#
# If `normalize_input=true` is provided, also normalizes the polynomials. 
# If `sort_input=true` is provided, also sorts the polynomials.
@timeit function f4_initialize_structs(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    normalize_input=true,
    sort_input=true
) where {M <: Monom, C <: Coeff}
    @log :debug "Initializing structs.."

    tablesize = hashtable_select_initial_size(ring, monoms)
    @log :debug "Initial hashtable size is $tablesize"

    # Basis for storing basis elements,
    # Pairset for storing critical pairs of basis elements,
    # Hashtable for hashing monomials stored in the basis
    basis = basis_initialize(ring, length(monoms), C)
    length(monoms) > typemax(Int16) && __basis_throw_maximum_size_exceeded(length(monoms))
    pairset = pairset_initialize(monom_entrytype(M))
    hashtable = hashtable_initialize(ring, params.rng, M, tablesize)

    # Filling the basis and hashtable with the given inputs
    basis_fill_data!(basis, hashtable, monoms, coeffs)
    hashtable_fill_divmasks!(hashtable)

    if params.changematrix
        basis_changematrix_initialize!(basis, hashtable)
    end

    if sort_input
        permutation = sort_polys_by_lead_increasing!(basis, hashtable, params.changematrix)
    else
        permutation = collect(1:(basis.nfilled))
    end

    # Divide each polynomial by the leading coefficient.
    # We do not need normalization in some cases, e.g., when computing the
    # normal forms
    if normalize_input
        basis_make_monic!(basis, params.arithmetic, params.changematrix)
    end

    basis, pairset, hashtable, permutation
end

# F4 reduction
@timeit function f4_reduction!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters
)
    # Re-enumerate matrix columns
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    # Call the linear algebra backend
    linalg_main!(matrix, basis, params)
    # Extract nonzero rows from the matrix into the basis
    matrix_convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht, params)
end

@timeit function f4_update!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M}
) where {M <: Monom}
    # total number of elements in the basis (old + new)
    npivs = basis.nfilled
    # number of potential critical pairs to add
    npairs = basis.nprocessed * npivs + div((npivs + 1) * npivs, 2)

    @invariant basis.nfilled >= basis.nprocessed
    @stat new_basis_elements = basis.nfilled - basis.nprocessed

    # make sure the pairset has enough space
    pairset_resize_if_needed!(pairset, npairs)
    pairset_size = length(pairset.pairs)

    # update pairset:
    # for each new element in basis..
    @inbounds for i in (basis.nprocessed + 1):(basis.nfilled)
        # ..check redundancy of new polynomial..
        basis_is_new_polynomial_redundant!(pairset, basis, ht, update_ht, i) && continue
        pairset_resize_lcms_if_needed!(pairset, basis.nfilled)
        # ..if not redundant, then add new S-pairs to pairset
        pairset_update!(pairset, basis, ht, update_ht, i)
    end

    basis_update!(basis, ht)
    pairset_size
end

@timeit function f4_symbolic_preprocessing!(
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    # The matrix already has rows that correspond to critical pairs. Here, we
    # find and add polynomial reducers to the matrix.
    #
    # Monomials that represent the columns of the matrix are stored in the
    # symbol_ht hashtable.
    symbol_load = symbol_ht.load
    ncols = matrix.ncols_left
    matrix_resize_upper_part_if_needed!(matrix, ncols + symbol_load)
    @log :debug "Finding reducers in the basis..." basis.nnonredundant

    # Traverse all monomials in symbol_ht and search for a polynomial reducer
    # for each monomial.
    #
    # Note that the size of hashtable grows as polynomials with new monomials
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
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hashval.hash, hashval.divmask, hashval.deg)
        matrix.ncols_left += 1
        f4_find_multiplied_reducer!(basis, matrix, ht, symbol_ht, MonomId(i))
        i += 1
    end

    # Shrink the matrix.
    resize!(matrix.upper_rows, matrix.nrows_filled_upper)

    nothing
end

# Performs autoreduction of basis elements inplace
function f4_autoreduce!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    params
) where {M}
    @log :debug "Entering autoreduction" basis

    etmp = monom_construct_const(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    matrix_reinitialize!(matrix, basis.nnonredundant)
    uprows = matrix.upper_rows

    # add all non redundant elements from the basis
    # as matrix upper rows
    @inbounds for i in 1:(basis.nnonredundant) #
        row_idx = matrix.nrows_filled_upper += 1
        uprows[row_idx] = matrix_polynomial_multiple_to_matrix_row!(
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
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)
    end

    # needed for correct column count in symbol hashtable
    matrix.ncols_left = matrix.ncols_left

    f4_symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        hv = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] = Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)
    end

    matrix_fill_column_to_monom_map!(matrix, symbol_ht)

    linalg_autoreduce!(matrix, basis, params)

    matrix_convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht, params)

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
end

function f4_select_tobereduced!(
    basis::Basis,
    tobereduced::Basis,
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable{M},
    ht::MonomialHashtable{M}
) where {M}

    # prepare to load all elems from tobereduced
    # to lower rows of the matrix
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
            matrix_polynomial_multiple_to_matrix_row!(matrix, symbol_ht, ht, h, etmp, gen)
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
    monomid::MonomId;
    sugar::Bool=false
)
    @inbounds monom = symbol_ht.monoms[monomid]
    @inbounds tmp = ht.monoms[1]
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

    # Compute the monomial quotient.
    success, tmp = monom_is_divisible!(tmp, monom, lead)
    if !success # division mask failed for some reason
        i += 1
        @goto Letsgo
    end
    # Using the fact that the hash is linear.
    @inbounds quotient_hash = symbol_ht.hashdata[monomid].hash - ht.hashdata[poly[1]].hash

    hashtable_resize_if_needed!(ht, length(poly))
    row = matrix_polynomial_multiple_to_matrix_row!(
        matrix,
        symbol_ht,
        ht,
        quotient_hash,
        tmp,
        poly,
        true
    )
    @inbounds row[1] = monomid
    @inbounds matrix.upper_rows[matrix.nrows_filled_upper + 1] = row
    @inbounds matrix.upper_to_coeffs[matrix.nrows_filled_upper + 1] = basis.nonredundant[i]

    # Insert the quotient into the main hashtable.
    # TODO: This line is here with one sole purpose -- to support tracing.
    #       Probably want to factor it out.
    @inbounds matrix.upper_to_mult[matrix.nrows_filled_upper + 1] =
        hashtable_insert!(ht, tmp)

    # Mark the current monomial as a pivot since a reducer has been found.
    @inbounds monom_hashval = symbol_ht.hashdata[monomid]
    @inbounds symbol_ht.hashdata[monomid] = Hashvalue(
        PIVOT_COLUMN,
        monom_hashval.hash,
        monom_hashval.divmask,
        monom_hashval.deg
    )
    matrix.nrows_filled_upper += 1

    nothing
end

# Returns N, the number of critical pairs of the smallest degree.
# Sorts the critical pairs so that the first N pairs are the smallest.
function pairset_lowest_degree_pairs!(pairset::Pairset)
    cnt = if true
        n_lowest_degree_pairs = pairset_partition_by_degree!(pairset)
        n_lowest_degree_pairs
    else
        sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
        pair_idx, _ = pairset_find_smallest_degree_pair(pairset)
        pair_idx
    end
    @invariant cnt > 0
    cnt
end

# Returns N, the number of critical pairs of the smallest sugar.
# Sorts the critical pairs so that the first N pairs are the smallest.
function pairset_lowest_sugar_pairs!(pairset::Pairset, sugar_cubes::Vector{SugarCube})
    sugar = sort_pairset_by_sugar!(pairset, 1, pairset.load - 1, sugar_cubes)
    @inbounds min_sugar = sugar[1]
    min_idx = 1
    @inbounds while min_idx < pairset.load && sugar[min_idx + 1] == min_sugar
        min_idx += 1
    end
    min_idx
end

# Discard all S-pairs of the lowest degree of lcm
# from the pairset
function f4_discard_normal!(
    pairset::Pairset,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable;
    maxpairs::Int=typemax(Int)
)
    npairs = pairset.load
    npairs = pairset_lowest_degree_pairs!(pairset)

    ps = pairset.pairs
    degs = pairset.degrees

    # if maxpairs is set or not
    if maxpairs != typemax(Int)
        sort_pairset_by_lcm!(pairset, npairs, ht)

        if npairs > maxpairs
            navailable = npairs
            npairs = maxpairs
            lastlcm = ps[npairs].lcm
            while npairs < navailable && ps[npairs + 1].lcm == lastlcm
                npairs += 1
            end
        end
    end

    @inbounds for i in 1:(pairset.load - npairs)
        ps[i] = ps[i + npairs]
        degs[i] = degs[i + npairs]
    end
    pairset.load -= npairs
end

# F4 critical pair selection.
function f4_select_critical_pairs!(
    pairset::Pairset{ExponentType},
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable;
    maxpairs::Int=typemax(Int),
    select_all::Bool=false
) where {ExponentType}
    # TODO: Why is this code type unstable !???
    npairs::Int = pairset.load
    if !select_all
        npairs = pairset_lowest_degree_pairs!(pairset)
    end
    npairs = min(npairs, maxpairs)
    @invariant npairs > 0
    ret[] += npairs

    ps = pairset.pairs
    degs = pairset.degrees
    @inbounds deg::ExponentType = degs[1]

    sort_pairset_by_lcm!(pairset, npairs, ht)

    # When `maxpairs` limits the number of selected pairs, we still add some
    # additional pairs which have the same lcm as the selected ones. 
    if npairs > maxpairs
        navailable = npairs
        npairs = maxpairs
        lastlcm = ps[npairs].lcm
        while npairs < navailable && ps[npairs + 1].lcm == lastlcm
            npairs += 1
        end
    end

    # At this stage, the first smallest pairs are selected. 
    # Add these pairs to the matrix.
    f4_add_critical_pairs_to_matrix!(pairset, npairs, basis, matrix, ht, symbol_ht)

    # Remove selected parirs from the pairset.
    @inbounds for i in 1:(pairset.load - npairs)
        ps[i] = ps[i + npairs]
        degs[i] = degs[i + npairs]
    end
    pairset.load -= npairs

    # @stat critical_pairs_deg = deg critical_pairs_count = npairs

    deg, npairs
end

# Adds the first `npairs` pairs from the pairset to the matrix.
# Write the corresponding monomials to symbolic hashtable.
#
# Assumes that the pairset is sorted w.r.t. the lcms of pairs.
function f4_add_critical_pairs_to_matrix!(
    pairset::Pairset,
    npairs::Int,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    matrix_reinitialize!(matrix, npairs)
    pairs = pairset.pairs
    uprows = matrix.upper_rows
    lowrows = matrix.lower_rows

    polys = Vector{Int}(undef, 2 * npairs)

    # For each pair...
    @inbounds etmp = ht.monoms[1]
    i = 1
    @inbounds while i <= npairs
        matrix.ncols_left += 1
        npolys = 1
        lcm = pairs[i].lcm
        j = i

        # ...collect all polynomials with the same lcm...
        while j <= npairs && pairs[j].lcm == lcm
            polys[npolys] = pairs[j].poly1
            npolys += 1
            polys[npolys] = pairs[j].poly2
            npolys += 1
            j += 1
        end
        npolys -= 1

        # ...and sort them by their position in the basis array.
        sort_generators_by_position!(polys, npolys)

        prev = polys[1]
        poly_monoms = basis.monoms[prev]
        vidx = poly_monoms[1]

        eidx = ht.monoms[vidx]
        elcm = ht.monoms[lcm]
        etmp = monom_division!(etmp, elcm, eidx)

        # hash of the quotient
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # The first polynomial from the pair is a reducer.
        row_idx = matrix.nrows_filled_upper += 1
        uprows[row_idx] = matrix_polynomial_multiple_to_matrix_row!(
            matrix,
            symbol_ht,
            ht,
            htmp,
            etmp,
            poly_monoms
        )
        matrix.upper_to_coeffs[row_idx] = prev
        matrix.upper_to_mult[row_idx] = hashtable_insert!(ht, etmp)

        # Mark the matrix column that corresponds to the lcm as a pivot.
        hv = symbol_ht.hashdata[uprows[row_idx][1]]
        symbol_ht.hashdata[uprows[row_idx][1]] =
            Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)

        # Add all polynomials with the same lcm to the lower part of matrix (to
        # be reduced).
        for k in 1:npolys
            # Duplicate polynomial. We can do so as long as the polynomials are
            # sorted.
            polys[k] == prev && continue

            # hashtable could have been reallocated, so refresh the pointer.
            elcm = ht.monoms[lcm]

            prev = polys[k]
            poly_monoms = basis.monoms[prev]
            vidx = poly_monoms[1]
            eidx = ht.monoms[vidx]
            etmp = monom_division!(etmp, elcm, eidx)
            # hash of the quotient
            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # The polynomial will be reduced.
            row_idx = matrix.nrows_filled_lower += 1
            lowrows[row_idx] = matrix_polynomial_multiple_to_matrix_row!(
                matrix,
                symbol_ht,
                ht,
                htmp,
                etmp,
                poly_monoms
            )
            matrix.lower_to_coeffs[row_idx] = prev
            matrix.lower_to_mult[row_idx] = hashtable_insert!(ht, etmp)

            hv = symbol_ht.hashdata[lowrows[row_idx][1]]
            symbol_ht.hashdata[lowrows[row_idx][1]] =
                Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)
        end

        i = j
    end

    resize!(matrix.lower_rows, matrix.nrows_filled_lower)
end

# move to basis.jl
function basis_well_formed(key, ring, basis, hashtable)
    if key in (:input_f4!, :input_f4_learn!, :input_f4_apply!)
        (isempty(basis.monoms) || isempty(basis.coeffs)) && return false
        (basis.size == 0 || basis.nfilled == 0) && return false
        !is_sorted_by_lead_increasing(basis, hashtable) && return false
    elseif key in (:output_f4!, :output_f4_learn!, :output_f4_apply!)
        !is_sorted_by_lead_increasing(basis, hashtable) && return false
        basis.nnonredundant ==
        length(basis.coeffs) ==
        length(basis.monoms) ==
        length(basis.divmasks) ==
        length(basis.nonredundant) ==
        length(basis.isredundant) || return false
        basis.nonredundant == collect(1:(basis.nnonredundant)) || return false
        any(!iszero, basis.isredundant) && return false
        any(c -> !isone(c[1]), basis.coeffs) && return false
    else
        return false
    end
    for i in 1:length(basis.coeffs)
        if !isassigned(basis.coeffs, i)
            if isassigned(basis.monoms, i)
                return false
            end
        else
            length(basis.coeffs[i]) == length(basis.monoms[i]) && continue
            if key in (:input_f4_apply!, :output_f4_apply!)
                @log :warn """
                key: $key
                Unlucky but perhaps not fatal cancellation in polynomial at index $(i) on apply stage.
                The number of monomials (expected): $(length(basis.monoms[i]))
                The number of monomials (got): $(length(basis.coeffs[i]))"""
            else
                return false
            end
        end
    end
    true
end

# F4 algorithm.
#
# Computes a groebner basis of the given `basis` inplace.
#
# Uses `pairset` to store critical pairs, 
# uses `hashtable` for hashing monomials,
# uses `tracer` to record information useful in subsequent runs.
#
# Input ivariants:
# - divmasks in the hashtable are set and correct,
# - basis is filled so that
#     basis.nfilled is the actual number of elements set,
#     basis.nprocessed     = 0,
#     basis.nnonredundant  = 0,
# - basis contains no zero polynomials (!!!).
#
# Output invariants:
# - basis.nprocessed == basis.nfilled == basis.nnonredundant
# - basis.monoms and basis.coeffs are of size basis.nprocessed
# - basis elements are sorted increasingly wrt the term ordering on lead elements
# - divmasks in basis are filled and coincide with divmasks in hashtable
@timeit function f4!(
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    tracer::TinyTraceF4,
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    # @invariant hashtable_well_formed(:input_f4!, ring, hashtable)
    @invariant basis_well_formed(:input_f4!, ring, basis, hashtable)
    # @invariant pairset_well_formed(:input_f4!, pairset, basis, ht)

    @log :debug "Entering F4."
    basis_make_monic!(basis, params.arithmetic, params.changematrix)

    matrix = matrix_initialize(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = hashtable_initialize_secondary(hashtable)
    symbol_ht = hashtable_initialize_secondary(hashtable)

    # add the first batch of critical pairs to the pairset
    @log :debug "Processing initial polynomials, generating first critical pairs"
    pairset_size = f4_update!(pairset, basis, hashtable, update_ht)
    update_tracer_pairset!(tracer, pairset_size)
    @log :debug "Out of $(basis.nfilled) polynomials, $(basis.nprocessed) are non-redundant"
    @log :debug "Generated $(pairset.load) critical pairs"

    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log :debug "F4: iteration $i"
        @log :debug "F4: available $(pairset.load) pairs"

        @log_memory_locals basis pairset hashtable update_ht symbol_ht

        # if the iteration is redundant according to the previous modular run
        if isready(tracer) && is_iteration_redundant(tracer, i)
            f4_discard_normal!(
                pairset,
                basis,
                matrix,
                hashtable,
                symbol_ht,
                maxpairs=params.maxpairs
            )
            matrix_reinitialize!(matrix, 0)
            hashtable_reinitialize!(symbol_ht)
            continue
        end

        # selects pairs for reduction from the pairset and puts these pairs into
        # the matrix rows
        f4_select_critical_pairs!(
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            maxpairs=params.maxpairs
        )
        @log :debug "After normal selection: available $(pairset.load) pairs"

        f4_symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)

        # reduces polys and obtains new potential basis elements
        f4_reduction!(ring, basis, matrix, hashtable, symbol_ht, params)

        update_tracer_iteration!(tracer, matrix.npivots == 0)

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = f4_update!(pairset, basis, hashtable, update_ht)
        update_tracer_pairset!(tracer, pairset_size)

        # clear symbolic hashtable
        # clear matrix
        matrix_reinitialize!(matrix, 0)
        hashtable_reinitialize!(symbol_ht)
        # matrix    = matrix_initialize(ring, C)
        # symbol_ht = hashtable_initialize_secondary(hashtable)

        if i > 10_000
            @log :warn "Something has gone wrong in F4. Error will follow."
            @log_memory_locals
            __f4_throw_maximum_iterations_exceeded(i)
        end
    end

    @stat f4_iterations = i

    set_ready!(tracer)
    set_final_basis!(tracer, basis.nfilled)

    if params.sweep
        @log :debug "Sweeping redundant elements in the basis"
        basis_sweep_redundant!(basis, hashtable)
    end

    basis_mark_redundant_elements!(basis)

    if params.reduced
        @log :debug "Autoreducing the final basis.."
        f4_autoreduce!(ring, basis, matrix, hashtable, symbol_ht, params)
    end

    basis_standardize!(
        ring,
        basis,
        hashtable,
        hashtable.ord,
        params.arithmetic,
        params.changematrix
    )

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    @invariant basis_well_formed(:output_f4!, ring, basis, hashtable)

    nothing
end

# Checks that all S-polynomials formed by the elements of the given basis reduce
# to zero.
@timeit function f4_isgroebner!(
    ring,
    basis::Basis{C},
    pairset,
    hashtable::MonomialHashtable{M},
    arithmetic::A
) where {M <: Monom, C <: Coeff, A <: AbstractArithmetic}
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(hashtable)
    update_ht = hashtable_initialize_secondary(hashtable)
    @log :debug "Forming S-polynomials"
    f4_update!(pairset, basis, hashtable, update_ht)
    isempty(pairset) && return true
    # Fill the F4 matrix
    f4_select_critical_pairs!(pairset, basis, matrix, hashtable, symbol_ht, select_all=true)
    f4_symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
    # Rename the columns and sort the rows of the matrix
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    # Reduce!
    linalg_isgroebner!(matrix, basis, arithmetic)
end

# Reduces each polynomial in the `tobereduced` by the polynomials from the `basis`.
@timeit function f4_normalform!(
    ring::PolyRing,
    basis::Basis{C},
    tobereduced::Basis{C},
    ht::MonomialHashtable,
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    matrix = matrix_initialize(ring, C)
    symbol_ht = hashtable_initialize_secondary(ht)
    # Fill the matrix
    f4_select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)
    f4_symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    matrix_fill_column_to_monom_map!(matrix, symbol_ht)
    # Reduce the matrix
    linalg_normalform!(matrix, basis, arithmetic)
    # Export the rows of the matrix back to the basis elements
    matrix_convert_rows_to_basis_elements_nf!(matrix, tobereduced, ht, symbol_ht)
    tobereduced
end
