# Main file that defines the f4! function.

# Functions in this file mostly accept a subset of these arguments:
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

# Performs gaussian row reduction of rows in the `matrix`
# and writes any nonzero results to `basis`
function reduction!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    linalg::Symbol,
    rng
)
    column_to_monom_mapping!(matrix, symbol_ht)

    sort_matrix_upper_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_lower_rows_increasing!(matrix) # for reduced, CD part

    linear_algebra!(ring, matrix, basis, linalg, rng)

    convert_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
end

function initialize_structs(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    normalize=true
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

    @log level = -4 "Sorting input polynomials by their leading terms in non-decreasing order"
    sort_polys_by_lead_increasing!(basis, hashtable)

    # Divide each polynomial by the leading coefficient
    if normalize
        @log level = -4 "Normalizing input polynomials"
        normalize_basis!(ring, basis)
    end

    basis, pairset, hashtable
end

function initialize_structs_learn(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters;
    normalize=true
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

    @log level = -4 "Sorting input polynomials by their leading terms in non-decreasing order"
    permutation = sort_polys_by_lead_increasing!(basis, hashtable)

    # Divide each polynomial by the leading coefficient
    if normalize
        @log level = -4 "Normalizing input polynomials"
        normalize_basis!(ring, basis)
    end

    @log level = -4 "Initializing computation graph"
    graph = initialize_computation_graph_f4(
        ring,
        deepcopy_basis(basis),
        basis,
        hashtable,
        permutation,
        params
    )

    graph, basis, pairset, hashtable
end

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
#
# Hashtable initial size is set to tablesize
function initialize_structs_no_normalize(
    ring::PolyRing,
    exponents::Vector{Vector{M}},
    coeffs_qq,
    coeffs_ff::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int
) where {M <: Monom, C <: Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)
    basis_ht = initialize_basis_hash_table(ring, rng, M, initial_size=tablesize)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs_ff)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # sort input, smaller leading terms first
    sort_polys_by_lead_increasing!(basis, basis_ht, coeffs_qq)

    basis, basis_ht
end

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
#
# Hashtable initial size is set to tablesize
function initialize_structs_ff(
    ring::PolyRing{Ch},
    exponents::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int
) where {Ch, M, C <: Coeff}
    coeffs_ff = [Vector{Ch}(undef, length(c)) for c in coeffs]
    initialize_structs_no_normalize(ring, exponents, coeffs, coeffs_ff, rng, tablesize)
end

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
function initialize_structs(
    ring::PolyRing,
    exponents::Vector{Vector{M}},
    coeffs_qq::Vector{Vector{T1}},
    coeffs_zz::Vector{Vector{T2}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int
) where {M, C <: Coeff, T1 <: CoeffQQ, T2 <: CoeffZZ}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)
    basis_ht = initialize_basis_hash_table(ring, rng, M, initial_size=tablesize)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # sort input, smaller leading terms first
    sort_polys_by_lead_increasing!(basis, basis_ht, coeffs_zz, coeffs_qq)

    # divide each polynomial by leading coefficient
    normalize_basis!(ring, basis)

    basis, basis_ht
end

# Initializes Basis with the given hashtable,
# fills input data from exponents and coeffs
function initialize_basis_using_existing_hashtable(
    ring::PolyRing,
    exponents::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable;
    sort_input=false,
    normalize_input=false
) where {M, C <: Coeff}
    basis = initialize_basis(ring, length(exponents), C)
    # fill the basis with the given inputs using the given hashtable as the
    # reference
    fill_data!(basis, present_ht, exponents, coeffs)
    if sort_input
        # sort input, smaller leading terms first
        @log level = -2 "Sorting input polynomials by the increasing leading term"
        sort_polys_by_lead_increasing!(basis, present_ht)
    end
    if normalize_input
        # divide each polynomial by leading coefficient
        # We do not need normalization for normal forms
        @log level = -2 "Normalizing input polynomials"
        normalize_basis!(ring, basis)
    end
    basis
end

# Initializes Basis with the given hashtable,
# fills input data from exponents and coeffs
function initialize_structs(
    ring::PolyRing,
    exponents::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int,
    present_ht::MonomialHashtable
) where {M, C <: Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, present_ht, exponents, coeffs)

    # sort input, smaller leading terms first
    sort_polys_by_lead_increasing!(basis, present_ht)

    # divide each polynomial by leading coefficient
    # We do not need normalization for normal forms
    # normalize_basis!(basis)

    basis, present_ht
end

# Initializes Basis with the given hashed exponents and coefficients
function initialize_structs(
    ring::PolyRing,
    hashedexps::Vector{Vector{MonomIdx}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable
) where {C <: Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, hashedexps, coeffs)
    basis.nfilled = length(hashedexps)

    # sort input, smaller leading terms first
    sort_polys_by_lead_increasing!(basis, present_ht)

    basis, present_ht
end

# Given a `basis` object that stores some groebner basis
# performs basis interreduction and writes the result to `basis` inplace
function reducegb_f4!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M}
) where {M}
    @log level = -5 "Entering autoreduction" basis

    etmp = construct_const_monom(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    reinitialize_matrix!(matrix, basis.nnonredundant)
    uprows = matrix.uprows

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

        matrix.up2coef[matrix.nrows] = basis.nonredundant[i]
        matrix.up2mult[matrix.nrows] = insert_in_hash_table!(ht, etmp)
        # set lead index as 1
        symbol_ht.hashdata[uprows[matrix.nrows][1]].idx = 1
    end

    # needed for correct column count in symbol hashtable
    matrix.ncols = matrix.nrows
    matrix.nup = matrix.nrows

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # set all pivots to unknown
    @inbounds for i in (symbol_ht.offset):(symbol_ht.load)
        symbol_ht.hashdata[i].idx = 1
    end

    column_to_monom_mapping!(matrix, symbol_ht)
    matrix.ncols = matrix.nleft + matrix.nright

    sort_matrix_upper_rows_decreasing!(matrix)

    exact_sparse_rref_interreduce!(ring, matrix, basis)

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
    resize!(matrix.lowrows, tobereduced.nfilled)

    etmp = construct_const_monom(M, ht.nvars)

    @inbounds for i in 1:(tobereduced.nfilled)
        matrix.nrows += 1
        gen = tobereduced.monoms[i]
        h = MonomHash(0)
        matrix.lowrows[matrix.nrows] =
            multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, gen)
        matrix.low2coef[matrix.nrows] = i
        matrix.low2mult[matrix.nrows] = insert_in_hash_table!(ht, etmp)
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
    vidx::MonomIdx
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

    # here found polynomial from basis with leading monom
    # dividing symbol_ht.monoms[vidx]
    if i <= basis.nnonredundant
        # reducers index and exponent in hash table
        @inbounds rpoly = basis.monoms[basis.nonredundant[i]]
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

        matrix.uprows[matrix.nup + 1] =
            multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, rpoly)
        @inbounds matrix.up2coef[matrix.nup + 1] = basis.nonredundant[i]
        # TODO: this line is here with the sole purpose -- to support tracing.
        # Probably want to factor it out.
        matrix.up2mult[matrix.nup + 1] = insert_in_hash_table!(ht, etmp)

        symbol_ht.hashdata[vidx].idx = 2
        matrix.nup += 1
        i += 1
    end

    nothing
end

# Recursively finds all polynomials from `basis` with the leading term
# that divides any of the monomials stored in hashtable `symbol_ht`,
# and writes all found polynomials to the `matrix`
function symbolic_preprocessing!(
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable
)
    symbol_load = symbol_ht.load

    nrr = matrix.ncols
    onrr = matrix.ncols

    while matrix.size <= nrr + symbol_load
        matrix.size *= 2
        resize!(matrix.uprows, matrix.size)
        resize!(matrix.up2coef, matrix.size)
        resize!(matrix.up2mult, matrix.size)
    end

    # for each lcm present in symbolic_ht set on select stage
    i = MonomIdx(symbol_ht.offset)
    #= First round, we add multiplied polynomials which divide =#
    #= a monomial exponent from selected spairs  =#
    @inbounds while i <= symbol_load
        # not a reducer
        if iszero(symbol_ht.hashdata[i].idx)
            symbol_ht.hashdata[i].idx = 1
            matrix.ncols += 1
            find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        end
        i += MonomIdx(1)
    end

    #= Second round, we add multiplied polynomials that divide  =#
    #= lcm added on previous for loop                            =#
    @inbounds while i <= symbol_ht.load
        if matrix.size == matrix.nup
            matrix.size *= 2
            resize!(matrix.uprows, matrix.size)
            resize!(matrix.up2coef, matrix.size)
            resize!(matrix.up2mult, matrix.size)
        end

        symbol_ht.hashdata[i].idx = 1
        matrix.ncols += 1
        find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        i += MonomIdx(1)
    end

    # shrink matrix sizes, set constants
    resize!(matrix.uprows, matrix.nup)

    matrix.nrows += matrix.nup - onrr
    matrix.nlow = matrix.nrows - matrix.nup
    matrix.size = matrix.nrows
end

# Returns the number of critical pairs of the smallest degree of lcm
function lowest_degree_pairs!(pairset::Pairset)
    sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
    ps = pairset.pairs
    @inbounds min_deg = ps[1].deg
    min_idx = 1
    @inbounds while min_idx < (pairset.load - 1) && ps[min_idx + 1].deg == min_deg
        min_idx += 1
    end
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

# Select all S-pairs of the lowest degree of lcm
# from the pairset and write the corresponding polynomials
# to the matrix
function select_normal!(
    pairset::Pairset,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable;
    maxpairs::Int=typemax(Int),
    selectall::Bool=false
)

    # number of selected pairs
    npairs = pairset.load
    if !selectall
        npairs = lowest_degree_pairs!(pairset)
    end
    ps = pairset.pairs

    npairs = min(npairs, maxpairs)
    # @info "Selected $(npairs) pairs"

    sort_pairset_by_lcm!(pairset, npairs, ht)

    if npairs > maxpairs
        navailable = npairs
        npairs = maxpairs
        lastlcm = ps[npairs].lcm
        while npairs < navailable && ps[npairs + 1].lcm == lastlcm
            npairs += 1
        end
    end

    @log level = -4 "Selected $(npairs) pairs in select_normal!"

    reinitialize_matrix!(matrix, npairs)

    uprows = matrix.uprows
    lowrows = matrix.lowrows

    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens = Vector{Int}(undef, 2 * npairs)

    deg = ps[1].deg

    # monomial buffer
    etmp = ht.monoms[1]
    i = 1
    @inbounds while i <= npairs
        matrix.ncols += 1
        load = 1
        lcm = ps[i].lcm
        j = i

        # we collect all generators with same lcm into gens
        while j <= npairs && ps[j].lcm == lcm
            gens[load] = ps[j].poly1
            load += 1
            gens[load] = ps[j].poly2
            load += 1
            j += 1
        end
        load -= 1

        # sort by the index in the basis (by=identity)
        sort_generators_by_position!(gens, load)

        # now we collect reducers and to-be-reduced polynomials

        # first generator index in groebner basis
        prev = gens[1]
        # first generator in hash table
        poly = basis.monoms[prev]
        # first generator lead monomial index in hash data
        vidx = poly[1]

        # first generator exponent
        eidx = ht.monoms[vidx]
        # exponent of lcm corresponding to first generator
        elcm = ht.monoms[lcm]
        etmp = monom_division!(etmp, elcm, eidx)
        # now etmp contents complement to eidx in elcm

        # hash of complement
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # add row as a reducer
        matrix.nup += 1
        uprows[matrix.nup] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
        # map upper row to index in basis
        matrix.up2coef[matrix.nup] = prev
        matrix.up2mult[matrix.nup] = insert_in_hash_table!(ht, etmp)

        # mark lcm column as reducer in symbolic hashtable
        symbol_ht.hashdata[uprows[matrix.nup][1]].idx = 2
        # increase number of rows set
        matrix.nrows += 1

        # over all polys with same lcm,
        # add them to the lower part of matrix
        @inbounds for k in 1:load
            # duplicate generator,
            # we can do so as long as generators are sorted
            if gens[k] == prev
                continue
            end

            # if the table was reallocated
            elcm = ht.monoms[lcm]

            # index in gb
            prev = gens[k]
            # poly of indices of monoms in hash table
            poly = basis.monoms[prev]
            vidx = poly[1]
            # leading monom idx
            eidx = ht.monoms[vidx]

            etmp = monom_division!(etmp, elcm, eidx)

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlow += 1
            lowrows[matrix.nlow] =
                multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
            # map lower row to index in basis
            matrix.low2coef[matrix.nlow] = prev
            matrix.low2mult[matrix.nlow] = insert_in_hash_table!(ht, etmp)

            symbol_ht.hashdata[lowrows[matrix.nlow][1]].idx = 2

            matrix.nrows += 1
        end

        i = j
    end

    resize!(matrix.lowrows, matrix.nrows - matrix.ncols)

    # remove selected parirs from pairset
    @inbounds for i in 1:(pairset.load - npairs)
        ps[i] = ps[i + npairs]
    end
    pairset.load -= npairs

    @log level = -3 "Selected $(npairs) pairs of degree $(deg) from pairset, $(pairset.load) pairs left"
    nothing
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
                @log level = 10^3 "Unlucky but probably not fatal cancellation at index $(i)" length(
                    basis.monoms[i]
                ) length(basis.coeffs[i])
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
    # TODO: here, decide on the number field arithmetic implementation
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

        @show_locals basis pairset hashtable update_ht symbol_ht

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
        select_normal!(
            pairset,
            basis,
            matrix,
            hashtable,
            symbol_ht,
            maxpairs=params.maxpairs
        )
        # Color with [F4]

        symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
        @log level = -100 "Formed a matrix of size X, DISPLAY_MATRIX"

        # reduces polys and obtains new potential basis elements
        reduction!(ring, basis, matrix, hashtable, symbol_ht, params.linalg, params.rng)
        @log level = -100 ""

        update_tracer_iteration!(tracer, matrix.npivots == 0)

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = update!(pairset, basis, hashtable, update_ht)
        update_tracer_pairset!(tracer, pairset_size)
        @log level = -100 "something something"

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
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht)
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
    hashtable::MonomialHashtable{M}
) where {M <: Monom, C <: Coeff}
    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hashtable(hashtable)
    update_ht = initialize_secondary_hashtable(hashtable)
    @log level = -3 "Forming S-polynomials"
    update!(pairset, basis, hashtable, update_ht)
    isempty(pairset) && return true
    # Fill the F4 matrix
    select_normal!(pairset, basis, matrix, hashtable, symbol_ht, selectall=true)
    symbolic_preprocessing!(basis, matrix, hashtable, symbol_ht)
    # Rename the columns and sort the rows of the matrix
    column_to_monom_mapping!(matrix, symbol_ht)
    sort_matrix_upper_rows_decreasing!(matrix)
    sort_matrix_lower_rows_increasing!(matrix)
    # Reduce!
    exact_sparse_rref_isgroebner!(ring, matrix, basis)
end

# Reduces each polynomial in the `tobereduced` by the polynomials from the `basis`.
function f4_normalform!(
    ring::PolyRing,
    basis::Basis{C},
    tobereduced::Basis{C},
    ht::MonomialHashtable
) where {C <: Coeff}
    matrix = initialize_matrix(ring, C)
    symbol_ht = initialize_secondary_hashtable(ht)
    # Fill the matrix
    select_tobereduced!(basis, tobereduced, matrix, symbol_ht, ht)
    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    column_to_monom_mapping!(matrix, symbol_ht)
    sort_matrix_upper_rows_decreasing!(matrix)
    # Reduce the matrix
    exact_sparse_rref_nf!(ring, matrix, tobereduced, basis)
    # Export the rows of the matrix back to the basis elements
    convert_rows_to_basis_elements_nf!(matrix, tobereduced, ht, symbol_ht)
    tobereduced
end
