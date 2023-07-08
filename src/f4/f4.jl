# Main file that defines the f4! function.

# Functions in this file mostly accept a subset of these arguments:
# ring - current polynomial ring,
# basis - a struct that stores polynomials,
# matrix - a struct that stores coefficients of polynomials to compute normal
#            forms,
# hashtable - a hashtable that stores monomials. (each monomial in the basis
#        points to a bucket in the hashtable)

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
    @log level = 3 "Initializing structs.."

    tablesize = select_hashtable_size(ring, monoms)
    @log level = 3 "Initial hashtable size is $tablesize"

    # Basis for storing basis elements,
    # Pairset for storing critical pairs of basis elements,
    # Hashtable for hashing monomials stored in the basis
    basis = initialize_basis(ring, length(monoms), C)
    pairset = initialize_pairset(powertype(M))
    hashtable = initialize_hashtable(ring, params.rng, M, tablesize)

    # Filling the basis and hashtable with the given inputs
    fill_data!(basis, hashtable, monoms, coeffs)
    fill_divmask!(hashtable)

    @log level = 4 "Hashtable:"

    @log level = 4 "Sorting input polynomials by their leading terms in non-decreasing order"
    sort_polys_by_lead_increasing!(basis, hashtable)

    # Divide each polynomial by the leading coefficient
    if normalize
        @log level = 4 "Normalizing input polynomials"
        normalize_basis!(ring, basis)
    end

    basis, pairset, hashtable
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
function initialize_structs_nf(
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
    # sort_polys_by_lead_increasing!(basis, present_ht)

    # divide each polynomial by leading coefficient
    # We do not need normalization for normal forms
    # normalize_basis!(basis)

    basis, present_ht
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
    basis.ntotal = length(hashedexps)

    # sort input, smaller leading terms first
    sort_polys_by_lead_increasing!(basis, present_ht)

    basis, present_ht
end

#------------------------------------------------------------------------------

# Given a `basis` object that stores some groebner basis
# performs basis interreduction and writes the result to `basis` inplace
function reducegb_f4!(
    ring::PolyRing,
    basis::Basis,
    matrix::MacaulayMatrix,
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M}
) where {M}
    etmp = make_zero_ev(M, ht.nvars)
    # etmp is now set to zero, and has zero hash

    reinitialize_matrix!(matrix, basis.ndivmasks)
    uprows = matrix.uprows

    # add all non redundant elements from basis
    # as matrix upper rows
    @inbounds for i in 1:(basis.ndivmasks) #
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

    basis.ntotal = matrix.npivots + basis.nprocessed
    basis.nprocessed = matrix.npivots

    # we may have added some multiples of reduced basis polynomials
    # from the matrix, so get rid of them
    k = 0
    i = 1
    @label Letsgo
    @inbounds while i <= basis.nprocessed
        @inbounds for j in 1:k
            if is_monom_divisible(
                basis.monoms[basis.ntotal - i + 1][1],
                basis.monoms[basis.nonredundant[j]][1],
                ht
            )
                i += 1
                @goto Letsgo
            end
        end
        k += 1
        basis.nonredundant[k] = basis.ntotal - i + 1
        basis.divmasks[k] = ht.hashdata[basis.monoms[basis.nonredundant[k]][1]].divmask
        i += 1
    end
    basis.ndivmasks = k
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
    reinitialize_matrix!(matrix, max(basis.ntotal, tobereduced.ntotal))
    resize!(matrix.lowrows, tobereduced.ntotal)

    etmp = make_zero_ev(M, ht.nvars)

    @inbounds for i in 1:(tobereduced.ntotal)
        matrix.nrows += 1
        gen = tobereduced.monoms[i]
        h = MonomHash(0)
        matrix.lowrows[matrix.nrows] =
            multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, gen)
        matrix.low2coef[matrix.nrows] = i
        matrix.low2mult[matrix.nrows] = insert_in_hash_table!(ht, etmp)
    end

    basis.ntotal
    basis.ndivmasks = basis.nprocessed = basis.ntotal
    basis.isredundant .= 0
    @inbounds for i in 1:(basis.ndivmasks)
        basis.nonredundant[i] = i
        basis.divmasks[i] = ht.hashdata[basis.monoms[i][1]].divmask
    end

    nothing
end

#------------------------------------------------------------------------------

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

    # searching for a poly from basis whose leading monom
    # divides the given exponent e
    i = 1
    @label Letsgo

    @inbounds while i <= basis.ndivmasks && (leaddiv[i] & ~divmask) != 0
        i += 1
    end

    # here found polynomial from basis with leading monom
    # dividing symbol_ht.monoms[vidx]
    if i <= basis.ndivmasks
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
        # TODO
        # TODO: insert in basis hashtable
        matrix.up2mult[matrix.nup + 1] = insert_in_hash_table!(ht, etmp)

        # up-size matrix
        symbol_ht.hashdata[vidx].idx = 2
        matrix.nup += 1
        i += 1
    end

    nothing
end

#------------------------------------------------------------------------------

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

    # TODO!
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
            # TODO::
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

#------------------------------------------------------------------------------

# Returns the number of critical pairs of the smallest degree of lcm
function lowest_degree_pairs!(pairset::Pairset)
    sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
    ps = pairset.pairs
    @inbounds min_deg = ps[1].deg
    min_idx = 0
    @inbounds while min_idx < pairset.load && ps[min_idx + 1].deg == min_deg
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
    maxpairs::Int=0
)
    npairs = lowest_degree_pairs!(pairset)
    # @debug "Discarded $(npairs) pairs"

    ps = pairset.pairs

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
        # TODO
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

    @log "Selected $(npairs) pairs of degree $(deg) from pairset, $(pairset.load) pairs left"
    nothing
end

function basis_well_formed(key, ring, basis, hashtable)
    if key in (:input_f4!, :input_f4_learn!, :input_f4_apply!)
        (isempty(basis.monoms) || isempty(basis.coeffs)) && return false
        (basis.size == 0 || basis.ntotal == 0) && return false
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
#     basis.ntotal is the actual number of set elements,
#     basis.nprocessed  = 0,
#     basis.ndivmasks  = 0,
# - basis contains no zero polynomials (!!!).
#
# Output invariants:
# - basis.nprocessed == basis.ntotal == basis.ndivmasks
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

    @log level = 5 "Entering F4."
    # TODO: decide on the number field arithmetic implementation
    normalize_basis!(ring, basis)
    
    matrix = initialize_matrix(ring, C)

    # initialize hash tables for update and symbolic preprocessing steps
    update_ht = initialize_secondary_hashtable(hashtable)
    symbol_ht = initialize_secondary_hashtable(hashtable)

    # makes basis fields valid,
    # does not copy,
    # checks for redundancy of new elems
    # TODO: move to basis.jl
    # plcm = Vector{MonomIdx}(undef, 0)
    # if isready(tracer)
    #     # TODO: `plcm` belongs into the pairset!!!
    #     resize!(plcm, final_basis_size(tracer) + 1)
    # end

    # add the first batch of critical pairs to the pairset
    @log level = 6 "Processing initial polynomials, generating first critical pairs"
    pairset_size = update!(pairset, basis, hashtable, update_ht)
    update_tracer_pairset!(tracer, pairset_size)
    @log level = 6 "Out of $(basis.ntotal) polynomials, $(basis.nprocessed) are non-redundant"
    @log level = 6 "Generated $(pairset.load) critical pairs"

    i = 0
    # While there are pairs to be reduced
    while !isempty(pairset)
        i += 1
        @log "F4: iteration $i"
        @log "F4: available $(pairset.load) pairs"

        # if the iteration is redundant according to the previous modular run
        # TODO: MOVE!
        if isready(tracer)
            if is_iteration_redundant(tracer, i)
                discard_normal!(pairset, basis, matrix, hashtable, symbol_ht)
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
        @log "Formed a matrix of size X, DISPLAY_MATRIX"

        # reduces polys and obtains new potential basis elements
        reduction!(ring, basis, matrix, hashtable, symbol_ht, params.linalg, params.rng)
        @log ""

        update_tracer_iteration!(tracer, matrix.npivots == 0)

        # update the current basis with polynomials produced from reduction,
        # does not copy,
        # checks for redundancy
        pairset_size = update!(pairset, basis, hashtable, update_ht)
        update_tracer_pairset!(tracer, pairset_size)
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

    set_ready!(tracer)
    set_final_basis!(tracer, basis.ntotal)

    # remove redundant elements
    filter_redundant!(basis)
    @log "Filtered elements marked redundant"

    if params.reduced
        @log "Autoreducing the final basis.."
        reducegb_f4!(ring, basis, matrix, hashtable, symbol_ht)
        @log "Autoreduced!"
    end

    standardize_basis!(ring, basis, hashtable, hashtable.ord)

    # @invariant hashtable_well_formed(:output_f4!, ring, hashtable)
    # @invariant basis_well_formed(:output_f4!, ring, basis, hashtable)

    nothing
end
