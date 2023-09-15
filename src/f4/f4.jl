# Main file that defines the f4! function.

# Functions here mostly accept or output a subset of these objects:
# ring - current polynomial ring,
# basis - a struct that stores polynomials,
# matrix - a struct that stores coefficients of polynomials to compute normal
#            forms,
# hashtable - a hashtable that stores monomials. (each monomial in the basis
#        points to a bucket in the hashtable)

@noinline __throw_maximum_iterations_exceeded(iters) =
    throw("""Something probably went wrong in Groebner.jl/F4. 
          The number of F4 iterations exceeded $iters. 
          Please consider submitting a GitHub issue.""")

# Given the polynomial ring and the arrays of monomials and coefficients,
# initializes and returns the following structures:
#   - basis: a Basis instance that stores polynomials,
#   - pairset: a Pairset instance that stores critical pairs,
#   - hashtable: a MonomialHashtable instance that stores monomials.
#   - permutation: a sorting permutation for input polynomials.
#
# If `normalize_input=true` is provided, normalizes the basis.
# If `sort_input=true` is provided, sorts the basis.
function initialize_structs(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    normalize_input=true,
    sort_input=true
) where {M <: Monom, C <: Coeff}
    @log level = -3 "Initializing structs.."

    tablesize = select_hashtable_size(ring, monoms)
    @log level = -3 "Initial hashtable size is $tablesize"

    # Basis for storing basis elements,
    # Pairset for storing critical pairs of basis elements,
    # Hashtable for hashing monomials stored in the basis
    basis = initialize_basis(ring, length(monoms), C)
    pairset = initialize_pairset(entrytype(M))
    hashtable = initialize_hashtable(ring, params.rng, M, tablesize)

    # Filling the basis and hashtable with the given inputs
    fill_data!(basis, hashtable, monoms, coeffs)
    fill_divmask!(hashtable)

    if sort_input
        permutation = sort_polys_by_lead_increasing!(basis, hashtable)
    else
        permutation = collect(1:(basis.nfilled))
    end

    # Divide each polynomial by the leading coefficient.
    # We do not need normalization in some cases, e.g., when computing the
    # normal forms
    if normalize_input
        normalize_basis!(ring, basis)
    end

    basis, pairset, hashtable, permutation
end

# Same as initialize_structs, but also initializes a computation graph
function initialize_structs_learn(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    normalize_input=true,
    sort_input=true
) where {M <: Monom, C <: Coeff}
    basis, pairset, hashtable, permutation = initialize_structs(
        ring,
        monoms,
        coeffs,
        params,
        normalize_input=normalize_input,
        sort_input=sort_input
    )

    graph = initialize_computation_graph_f4(
        ring,
        deepcopy_basis(basis),
        basis,
        hashtable,
        permutation,
        params
    )

    graph, basis, pairset, hashtable, permutation
end

# Same as initialize_structs, but uses an existing hashtable
function initialize_basis_using_existing_hashtable(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable;
) where {M, C <: Coeff}
    basis = initialize_basis(ring, length(monoms), C)
    fill_data!(basis, present_ht, monoms, coeffs)
    basis
end

# F4 reduction
@timed_block function reduction!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters
)
    # Construct a mapping from monomials to matrix columns and re-enumerate
    # matrix columns
    column_to_monom_mapping!(matrix, symbol_ht)
    linear_algebra!(matrix, basis, params)
    # Extract nonzero rows from the matrix into the basis
    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
end

# F4 symbolic preprocessing
@timed_block function symbolic_preprocessing!(
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    # 1. The matrix already has rows added on the critical pair selection stage.
    #    Here, we find and add polynomial reducers to the matrix.
    # 2. Monomials that represent the columns of the matrix are stored in the
    #    symbol_ht hashtable.
    symbol_load = symbol_ht.load
    ncols = matrix.ncolumns

    resize_matrix_upper_part_if_needed!(matrix, ncols + symbol_load)

    # 3. Traverse all monomials in symbol_ht and search for a polynomial reducer
    #    for each monomial.
    # NOTE: note that the size of hashtable grows as polynomials with new
    # monomials are added to the matrix, and the loop accounts for that
    i = MonomIdx(symbol_ht.offset)
    @inbounds while i <= symbol_ht.load
        if symbol_ht.hashdata[i].idx != NON_PIVOT_COLUMN
            i += MonomIdx(1)
            continue
        end
        resize_matrix_upper_part_if_needed!(matrix, matrix.nupper + 1)

        hashval = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] =
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hashval.hash, hashval.divmask, hashval.deg)
        matrix.ncolumns += 1
        find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        i += MonomIdx(1)
    end

    # Shrink the matrix, set the dimensions
    resize!(matrix.upper_rows, matrix.nupper)
    matrix.nrows += matrix.nupper - ncols
    matrix.nlower = matrix.nrows - matrix.nupper
    matrix.size = matrix.nrows
end

# Given a `basis` object that stores some groebner basis
# performs basis interreduction and writes the result to `basis` inplace
function reducegb_f4!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    params
) where {M}
    @log level = -7 "Entering autoreduction" basis

    etmp = construct_const_monom(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    reinitialize_matrix!(matrix, basis.nnonredundant)
    uprows = matrix.upper_rows

    # add all non redundant elements from the basis
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
        hv = symbol_ht.hashdata[uprows[matrix.nrows][1]]
        symbol_ht.hashdata[uprows[matrix.nrows][1]] =
            Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)
    end

    # needed for correct column count in symbol hashtable
    matrix.ncolumns = matrix.nrows
    matrix.nupper = matrix.nrows

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        hv = symbol_ht.hashdata[i]
        symbol_ht.hashdata[i] = Hashvalue(UNKNOWN_PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)
    end

    column_to_monom_mapping!(matrix, symbol_ht)
    matrix.ncolumns = matrix.nleft + matrix.nright

    linear_algebra_reducegb!(matrix, basis, params)

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)

    basis.nfilled = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

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
end

function select_tobereduced!(
    basis::Basis,
    tobereduced::Basis,
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable{M},
    ht::MonomialHashtable{M}
) where {M}

    # prepare to load all elems from tobereduced
    # to lower rows of the matrix
    reinitialize_matrix!(matrix, max(basis.nfilled, tobereduced.nfilled))
    resize!(matrix.lower_rows, tobereduced.nfilled)
    resize!(matrix.coeffs, tobereduced.nfilled)

    etmp = construct_const_monom(M, ht.nvars)

    @inbounds for i in 1:(tobereduced.nfilled)
        matrix.nrows += 1
        gen = tobereduced.monoms[i]
        h = MonomHash(0)
        matrix.lower_rows[matrix.nrows] =
            multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, gen)
        matrix.lower_to_coeffs[matrix.nrows] = i
        # TODO: not really needed here
        matrix.lower_to_mult[matrix.nrows] = insert_in_hash_table!(ht, etmp)
        matrix.coeffs[matrix.nrows] = tobereduced.coeffs[i]
    end

    basis.nfilled
    basis.nnonredundant = basis.nprocessed = basis.nfilled
    basis.isredundant .= 0
    @inbounds for i in 1:(basis.nnonredundant)
        basis.nonredundant[i] = i
        basis.divmasks[i] = ht.hashdata[basis.monoms[i][1]].divmask
    end

    nothing
end

# Finds a polynomial from the `basis` 
# with leading term that divides monomial `vidx`. 
# If such polynomial was found, 
# writes the divisor polynomial to the hashtable `symbol_ht`
function find_multiplied_reducer!(
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    vidx::MonomIdx;
    sugar::Bool=false
)
    e = symbol_ht.monoms[vidx]
    etmp = ht.monoms[1]
    divmask = symbol_ht.hashdata[vidx].divmask
    leaddiv = basis.divmasks

    # Searching for a poly from basis whose leading monom divides the given
    # exponent e
    i = 1
    @label Letsgo

    @inbounds while i <= basis.nnonredundant
        # TODO: rethink division masks to support more variables
        if ht.use_divmask && is_divmask_divisible(divmask, leaddiv[i])
            break
        else
            e2 = ht.monoms[basis.monoms[basis.nonredundant[i]][1]]
            if is_monom_divisible(e, e2)
                break
            end
        end
        i += 1
    end

    # Reducer is not found, yield
    i > basis.nnonredundant && return nothing

    # Here, found polynomial from basis with leading monom
    # dividing symbol_ht.monoms[vidx]

    # reducers index and exponent in hash table
    @inbounds rpoly = basis.monoms[basis.nonredundant[i]]
    resize_hashtable_if_needed!(ht, length(rpoly))

    @inbounds rexp = ht.monoms[rpoly[1]]

    # precisely, etmp = e .- rexp 
    flag, etmp = is_monom_divisible!(etmp, e, rexp)
    if !flag
        i += 1
        @goto Letsgo
    end
    # now etmp = e // rexp in terms of monomias,
    # (!) hash is linear
    @inbounds h = symbol_ht.hashdata[vidx].hash - ht.hashdata[rpoly[1]].hash

    matrix.upper_rows[matrix.nupper + 1] =
        multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, rpoly)
    @inbounds matrix.upper_to_coeffs[matrix.nupper + 1] = basis.nonredundant[i]
    # TODO: this line is here with the sole purpose -- to support tracing.
    # Probably want to factor it out.
    matrix.upper_to_mult[matrix.nupper + 1] = insert_in_hash_table!(ht, etmp)
    # if sugar
    #     # updates sugar
    #     poly = basis.nonredundant[i]
    #     new_poly_sugar = totaldeg(etmp) + basis.sugar_cubes[poly]
    #     matrix.upper_to_sugar[matrix.nupper + 1] = new_poly_sugar
    # end

    hv = symbol_ht.hashdata[vidx]
    symbol_ht.hashdata[vidx] = Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)

    matrix.nupper += 1
    i += 1

    nothing
end

# Returns N, the number of critical pairs of the smallest degree.
# Sorts the critical pairs so that the first N pairs are the smallest.
function lowest_degree_pairs!(pairset::Pairset)
    sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
    ps = pairset.pairs
    @inbounds min_deg = ps[1].deg
    min_idx = 1
    @inbounds while min_idx < pairset.load && ps[min_idx + 1].deg == min_deg
        min_idx += 1
    end
    min_idx
end

# Returns N, the number of critical pairs of the smallest sugar.
# Sorts the critical pairs so that the first N pairs are the smallest.
function lowest_sugar_pairs!(pairset::Pairset, sugar_cubes::Vector{SugarCube})
    @log level = -1 "Sugar cubes" sugar_cubes
    sugar = sort_pairset_by_sugar!(pairset, 1, pairset.load - 1, sugar_cubes)
    @inbounds min_sugar = sugar[1]
    min_idx = 1
    @inbounds while min_idx < pairset.load && sugar[min_idx + 1] == min_sugar
        min_idx += 1
    end
    @log level = -1 "Selected pairs sugar" sugar min_idx min_sugar
    min_idx
end

# Discard all S-pairs of the lowest degree of lcm
# from the pairset
function discard_normal!(
    pairset::Pairset,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable;
    maxpairs::Int=typemax(Int)
)
    npairs = pairset.load
    npairs = lowest_degree_pairs!(pairset)
    # @debug "Discarded $(npairs) pairs"

    ps = pairset.pairs

    # if maxpairs is set
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

    @debug "Discarded $(npairs) pairs"

    @inbounds for i in 1:(pairset.load - npairs)
        ps[i] = ps[i + npairs]
    end
    pairset.load -= npairs
end

# F4 critical pair selection.
#
# Selects critical pairs from the pairset and adds the corresponding polynomials
# to the matrix
#
# If `maxpairs=N` is provided, the number of critical pairs is limited by `N`
# (modulo some technical details).
# If `select_all=true` is provided, selects all critical pairs.
@timed_block function select_critical_pairs!(
    pairset::Pairset,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    selection_strategy::Symbol;
    maxpairs::Int=typemax(Int),
    select_all::Bool=false
)
    # Here, the following happens.
    # 1. The pairset is sorted according to the given selection strategy and the
    #    number of selected critical pairs is decided.

    npairs = pairset.load
    if !select_all
        if selection_strategy === :normal
            npairs = lowest_degree_pairs!(pairset)
        else
            @assert selection_strategy === :sugar
            npairs = lowest_sugar_pairs!(pairset, basis.sugar_cubes)
        end
    end
    npairs = min(npairs, maxpairs)
    @assert npairs > 0
    ps = pairset.pairs
    deg = ps[1].deg

    # 2. Selected pairs in the pairset are sorted once again, now with respect
    #    to a monomial ordering on the LCMs
    sort_pairset_by_lcm!(pairset, npairs, ht)

    # NOTE: when `maxpairs` limits the number of selected pairs, we still add
    # some additional pairs which have the same lcm as the selected ones 
    if npairs > maxpairs
        navailable = npairs
        npairs = maxpairs
        lastlcm = ps[npairs].lcm
        while npairs < navailable && ps[npairs + 1].lcm == lastlcm
            npairs += 1
        end
    end

    # 3. At this stage, we know that the first `npairs` pairs in the pairset are 
    #    selected. We add these pairs to the matrix
    @log level = -4 "Selected $(npairs) critical pairs"
    add_critical_pairs_to_matrix!(pairset, npairs, basis, matrix, ht, symbol_ht)

    # 4. Remove selected parirs from the pairset
    @inbounds for i in 1:(pairset.load - npairs)
        ps[i] = ps[i + npairs]
    end
    pairset.load -= npairs

    @log level = -3 "Selected $(npairs) pairs of degree $(deg) from pairset, $(pairset.load) pairs left"
    nothing
end

# Adds the first `npairs` pairs from the pairset to the matrix
function add_critical_pairs_to_matrix!(
    pairset::Pairset,
    npairs::Int,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    # 
    #
    #
    #
    reinitialize_matrix!(matrix, npairs)
    pairs = pairset.pairs
    uprows = matrix.upper_rows
    lowrows = matrix.lower_rows

    polys = Vector{Int}(undef, 2 * npairs)
    # monomial buffer
    etmp = ht.monoms[1]
    i = 1
    @inbounds while i <= npairs
        matrix.ncolumns += 1
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

        # sort by the index in the basis (by=identity)
        sort_generators_by_position!(polys, npolys)

        # now we collect reducers and to-be-reduced polynomials

        # first generator index in groebner basis
        prev = polys[1]
        # first generator in hash table
        poly_monoms = basis.monoms[prev]
        # first generator lead monomial index in hash data
        vidx = poly_monoms[1]

        # first generator exponent
        eidx = ht.monoms[vidx]
        # exponent of lcm corresponding to first generator
        elcm = ht.monoms[lcm]
        etmp = monom_division!(etmp, elcm, eidx)
        # now etmp contents complement to eidx in elcm

        # hash of complement
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # add row as a reducer
        matrix.nupper += 1
        uprows[matrix.nupper] =
            multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly_monoms)
        # map upper row to index in basis
        matrix.upper_to_coeffs[matrix.nupper] = prev
        matrix.upper_to_mult[matrix.nupper] = insert_in_hash_table!(ht, etmp)

        # mark lcm column as reducer in symbolic hashtable
        hv = symbol_ht.hashdata[uprows[matrix.nupper][1]]
        symbol_ht.hashdata[uprows[matrix.nupper][1]] =
            Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)

        # increase number of rows set
        matrix.nrows += 1

        # over all polys with same lcm,
        # add them to the lower part of matrix
        for k in 1:npolys
            # duplicate generator,
            # we can do so as long as generators are sorted
            if polys[k] == prev
                continue
            end

            # if the table was reallocated
            elcm = ht.monoms[lcm]

            # index in gb
            prev = polys[k]
            # poly of indices of monoms in hash table
            poly_monoms = basis.monoms[prev]
            vidx = poly_monoms[1]
            # leading monom idx
            eidx = ht.monoms[vidx]

            etmp = monom_division!(etmp, elcm, eidx)

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlower += 1
            lowrows[matrix.nlower] =
                multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly_monoms)
            # map lower row to index in basis
            matrix.lower_to_coeffs[matrix.nlower] = prev
            matrix.lower_to_mult[matrix.nlower] = insert_in_hash_table!(ht, etmp)

            hv = symbol_ht.hashdata[lowrows[matrix.nlower][1]]
            symbol_ht.hashdata[lowrows[matrix.nlower][1]] =
                Hashvalue(PIVOT_COLUMN, hv.hash, hv.divmask, hv.deg)

            matrix.nrows += 1
        end

        i = j
    end

    resize!(matrix.lower_rows, matrix.nrows - matrix.ncolumns)
end

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
                @log level = 10^3 """
                Unlucky but probably not fatal cancellation in polynomial at index $(i) on apply stage.
                The number of monomials: $(length(basis.monoms[i]))
                The number of coefficients: $(length(basis.coeffs[i]))"""
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
# - divmasks in the ht are set,
# - basis is filled so that
#     basis.nfilled is the actual number of set elements,
#     basis.nprocessed  = 0,
#     basis.nnonredundant  = 0,
# - basis contains no zero polynomials (!!!).
#
# Output invariants:
# - basis.nprocessed == basis.nfilled == basis.nnonredundant
# - basis.monoms and basis.coeffs are of size basis.nprocessed
# - basis elements are sorted increasingly wrt the term ordering on lead elements
# - divmasks in basis are filled and coincide with divmasks in hashtable
function f4!(
    ring::PolyRing,
    basis::Basis{C},
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    tracer::Tracer,
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    # @invariant hashtable_well_formed(:input_f4!, ring, hashtable)
    @invariant basis_well_formed(:input_f4!, ring, basis, hashtable)
    # @invariant pairset_well_formed(:input_f4!, pairset, basis, ht)

    @log level = -2 "Entering F4."
    normalize_basis!(ring, basis)

    matrix = initialize_matrix(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = initialize_secondary_hashtable(hashtable)
    symbol_ht = initialize_secondary_hashtable(hashtable)

    # add the first batch of critical pairs to the pairset
    @log level = -3 "Processing initial polynomials, generating first critical pairs"
    pairset_size = update!(pairset, basis, hashtable, update_ht)
    update_tracer_pairset!(tracer, pairset_size)
    @log level = -3 "Out of $(basis.nfilled) polynomials, $(basis.nprocessed) are non-redundant"
    @log level = -3 "Generated $(pairset.load) critical pairs"

    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log level = -3 "F4: iteration $i"
        @log level = -3 "F4: available $(pairset.load) pairs"

        @log_memory_locals basis pairset hashtable update_ht symbol_ht

        # if the iteration is redundant according to the previous modular run
        if isready(tracer)
            if is_iteration_redundant(tracer, i)
                discard_normal!(
                    pairset,
                    basis,
                    matrix,
                    hashtable,
                    symbol_ht,
                    maxpairs=params.maxpairs
                )
                matrix    = initialize_matrix(ring, C)
                symbol_ht = initialize_secondary_hashtable(hashtable)
                continue
            end
        end

        # selects pairs for reduction from pairset following normal strategy
        # (minimal lcm degrees are selected),
        # and puts these into the matrix rows
        select_critical_pairs!(
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            params.selection_strategy,
            maxpairs=params.maxpairs
        )
        @log level = -3 "After normal selection: available $(pairset.load) pairs"
        @log level = -3 repr_basis(basis)

        symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)

        # reduces polys and obtains new potential basis elements
        reduction!(ring, basis, matrix, hashtable, symbol_ht, params)

        update_tracer_iteration!(tracer, matrix.npivots == 0)

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = update!(pairset, basis, hashtable, update_ht)
        update_tracer_pairset!(tracer, pairset_size)

        # clear symbolic hashtable
        # clear matrix
        matrix    = initialize_matrix(ring, C)
        symbol_ht = initialize_secondary_hashtable(hashtable)

        if i > 10_000
            @log level = 1 "Something has gone wrong in F4. Error will follow."
            @log_memory_locals
            __throw_maximum_iterations_exceeded(i)
        end
    end

    set_ready!(tracer)
    set_final_basis!(tracer, basis.nfilled)

    if params.sweep
        @log level = -3 "Sweeping redundant elements in the basis"
        sweep_redundant!(basis, hashtable)
    end

    # mark redundant elements
    mark_redundant!(basis)
    @log level = -3 "Filtered elements marked redundant"

    if params.reduced
        @log level = -2 "Autoreducing the final basis.."
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht, params)
        @log level = -3 "Autoreduced!"
    end

    standardize_basis!(ring, basis, hashtable, hashtable.ord)

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    @invariant basis_well_formed(:output_f4!, ring, basis, hashtable)

    nothing
end

# Checks that all S-polynomials formed by the elements of the given basis reduce
# to zero.
function f4_isgroebner!(
    ring,
    basis::Basis{C},
    pairset,
    hashtable::MonomialHashtable{M},
    arithmetic::A
) where {M <: Monom, C <: Coeff, A <: AbstractArithmetic}
    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hashtable(hashtable)
    update_ht = initialize_secondary_hashtable(hashtable)
    @log level = -3 "Forming S-polynomials"
    update!(pairset, basis, hashtable, update_ht)
    isempty(pairset) && return true
    # Fill the F4 matrix
    select_critical_pairs!(
        pairset,
        basis,
        matrix,
        hashtable,
        symbol_ht,
        :normal,
        select_all=true
    )
    symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
    # Rename the columns and sort the rows of the matrix
    column_to_monom_mapping!(matrix, symbol_ht)
    # Reduce!
    linear_algebra_isgroebner!(matrix, basis, arithmetic)
end

# Reduces each polynomial in the `tobereduced` by the polynomials from the `basis`.
function f4_normalform!(
    ring::PolyRing,
    basis::Basis{C},
    tobereduced::Basis{C},
    ht::MonomialHashtable,
    arithmetic::A
) where {C <: Coeff, A <: AbstractArithmetic}
    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hashtable(ht)
    # Fill the matrix
    select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)
    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    column_to_monom_mapping!(matrix, symbol_ht)
    # Reduce the matrix
    linear_algebra_normalform!(matrix, basis, arithmetic)
    # Export the rows of the matrix back to the basis elements
    convert_rows_to_basis_elements_nf!(matrix, tobereduced, ht, symbol_ht)
    tobereduced
end
