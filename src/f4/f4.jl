
#=
    The file contains implementation of the F4 algorithm,
    described originally by Jean-Charles Faugere in
        "A new efficient algorithm for computing Grobner bases"
=#

#------------------------------------------------------------------------------

# Prepares the given set of polynomials L to be reduced by the ideal G,
# starts the reduction routine,
# and returns reduced polynomials
function reduction!(
    basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable,
    linalg::Symbol,
    rng)

    convert_hashes_to_columns!(matrix, symbol_ht)

    sort_matrix_rows_decreasing!(matrix) # for pivots,  AB part
    sort_matrix_rows_increasing!(matrix) # for reduced, CD part

    # exact_sparse_linear_algebra!(matrix, basis)
    linear_algebra!(matrix, basis, Val(linalg), rng)

    convert_matrix_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
end

#------------------------------------------------------------------------------

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
#
# Hashtable initial size is set to tablesize
function initialize_structures(
    ring::PolyRing,
    exponents::Vector{Vector{ExponentVector}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int) where {C<:Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)
    basis_ht = initialize_basis_hash_table(ring, rng, initial_size=tablesize)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # sort input, smaller leading terms first
    sort_gens_by_lead_increasing!(basis, basis_ht)

    # divide each polynomial by leading coefficient
    normalize_basis!(basis)

    basis, basis_ht
end

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
#
# Hashtable initial size is set to tablesize
function initialize_structures_no_normalize(
    ring::PolyRing,
    exponents::Vector{Vector{ExponentVector}},
    coeffs_qq,
    coeffs_ff::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int) where {C<:Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)
    basis_ht = initialize_basis_hash_table(ring, rng, initial_size=tablesize)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs_ff)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # sort input, smaller leading terms first
    sort_gens_by_lead_increasing!(basis, basis_ht, coeffs_qq)

    basis, basis_ht
end

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
#
# Hashtable initial size is set to tablesize
function initialize_structures_ff(
    ring::PolyRing,
    exponents::Vector{Vector{ExponentVector}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int) where {C<:Coeff}

    coeffs_ff = [Vector{UInt64}(undef, length(c)) for c in coeffs]
    initialize_structures_no_normalize(ring, exponents, coeffs, coeffs_ff, rng, tablesize)
end

# Initializes Basis and MonomialHashtable structures,
# fills input data from exponents and coeffs
function initialize_structures(
    ring::PolyRing,
    exponents::Vector{Vector{ExponentVector}},
    coeffs_qq::Vector{Vector{T1}},
    coeffs_zz::Vector{Vector{T2}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int) where {C<:Coeff,T1<:CoeffQQ,T2<:CoeffZZ}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)
    basis_ht = initialize_basis_hash_table(ring, rng, initial_size=tablesize)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, basis_ht, exponents, coeffs)

    # every monomial in hashtable is associated with its divmask
    # to perform divisions faster. Filling those
    fill_divmask!(basis_ht)

    # sort input, smaller leading terms first
    sort_gens_by_lead_increasing!(basis, basis_ht, coeffs_zz, coeffs_qq)

    # divide each polynomial by leading coefficient
    normalize_basis!(basis)

    basis, basis_ht
end


# Initializes Basis with the given hashtable,
# fills input data from exponents and coeffs
function initialize_structures(
    ring::PolyRing,
    exponents::Vector{Vector{ExponentVector}},
    coeffs::Vector{Vector{C}},
    rng::Random.AbstractRNG,
    tablesize::Int,
    present_ht::MonomialHashtable) where {C<:Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, length(exponents), C)

    # filling the basis and hashtable with the given inputs
    fill_data!(basis, present_ht, exponents, coeffs)

    # sort input, smaller leading terms first
    sort_gens_by_lead_increasing!(basis, present_ht)

    # divide each polynomial by leading coefficient
    # We do not need normalization for normal forms
    # normalize_basis!(basis)

    basis, present_ht
end

# Initializes Basis with the given hashed exponents and coefficients
function initialize_structures(
    ring::PolyRing,
    hashedexps::Vector{Vector{ExponentIdx}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable) where {C<:Coeff}

    # basis for storing basis elements,
    # pairset for storing critical pairs of basis elements to assess,
    # hashtable for hashing monomials occuring in the basis
    basis = initialize_basis(ring, hashedexps, coeffs)
    basis.ntotal = length(hashedexps)

    # sort input, smaller leading terms first
    sort_gens_by_lead_increasing!(basis, present_ht)

    basis, present_ht
end

#------------------------------------------------------------------------------

function reducegb_relaxed_f4!(
    basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    for i in 1:basis.nlead
        etmp = ht.exponents[1]
        etmp = zero(etmp)
        # etmp is now set to zero, and has a zero hash

        reinitialize_matrix!(matrix, basis.nlead)
        uprows = matrix.uprows

        #=
        @warn "entering reduce"
        dump(basis, maxdepth=5)
        @warn "ht"
        println(ht.exponents[1:10])
        @warn "symbol"
        println(symbol_ht.exponents[1:10])
        println("NONRED ", basis.nonred)
        =#

        # add all non redundant elements from basis
        # as matrix upper rows

        matrix.nrows += 1
        uprows[matrix.nrows] = multiplied_poly_to_matrix_row!(
            symbol_ht, ht, UInt32(0), etmp,
            basis.gens[basis.nonred[i]])

        matrix.up2coef[matrix.nrows] = basis.nonred[i]
        # set lead index as 1
        symbol_ht.hashdata[uprows[matrix.nrows][1]].idx = 1

        # needed for correct counting in symbol
        matrix.ncols = matrix.nrows
        matrix.nup = matrix.nrows

        #=
        @warn "after multiplied_poly_to_matrix_row"
        dump(basis, maxdepth=5)
        @warn "matrix"
        dump(matrix, maxdepth=5)
        @warn "ht"
        println(ht.exponents[1:10])
        =#

        symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
        # all pivots are unknown
        for j in symbol_ht.offset:symbol_ht.load
            symbol_ht.hashdata[j].idx = 1
        end

        #=
        @warn "after symbolic"
        dump(basis, maxdepth=5)
        @warn "matrix"
        dump(matrix, maxdepth=5)
        @warn "ht"
        println(ht.exponents[1:10])
        =#

        # x1*x2^6 + 413612941*x1*x2^3
        #  x1^2*x2^3 + 174101409*x1*x2^5
        convert_hashes_to_columns!(matrix, symbol_ht)
        matrix.ncols = matrix.nleft + matrix.nright

        #=
        @warn "after convert"
        dump(matrix, maxdepth=5)
        =#

        sort_matrix_rows_decreasing!(matrix)

        println(matrix.ncols)

        println(basis.ndone)

        interreduce_matrix_rows!(matrix, basis)

        #=
        @warn "after interreduce"
        dump(matrix, maxdepth=5)
        =#
        # TODO
        # convert_matrix_rows_to_basis_elements_use_symbol!(matrix, basis)
        convert_matrix_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
        # no longer need in two hashtables
        # TODO
        # ht = symbol_ht

        basis.ntotal = matrix.npivots + basis.ndone
        basis.ndone = matrix.npivots

        #=
        @warn "basis"
        dump(basis, maxdepth=5)
        =#

        #= we may have added some multiples of reduced basis polynomials
        * from the matrix, so we get rid of them. =#
        k = 0
        i = 1
        @label Letsgo
        while i <= basis.ndone
            @inbounds for j in 1:k
                if is_monom_divisible(
                    basis.gens[basis.ntotal-i+1][1],
                    basis.gens[basis.nonred[j]][1],
                    ht)

                    i += 1
                    @goto Letsgo
                end
            end
            k += 1
            basis.nonred[k] = basis.ntotal - i + 1
            # xd
            basis.lead[k] = ht.hashdata[basis.gens[basis.nonred[k]][1]].divmask
            i += 1
        end
        basis.nlead = k

    end

    # TODO
    # sort_gens_by_lead_increasing_in_reduce!(basis, ht)
end

function reducegb_f4!(
    basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    etmp = ht.exponents[1]
    etmp = zero(etmp)
    # etmp is now set to zero, and has a zero hash

    reinitialize_matrix!(matrix, basis.nlead)
    uprows = matrix.uprows

    #=
    @warn "entering reduce"
    dump(basis, maxdepth=5)
    @warn "ht"
    println(ht.exponents[1:10])
    @warn "symbol"
    println(symbol_ht.exponents[1:10])
    println("NONRED ", basis.nonred)
    =#

    # add all non redundant elements from basis
    # as matrix upper rows
    @inbounds for i in 1:basis.nlead #
        matrix.nrows += 1
        uprows[matrix.nrows] = multiplied_poly_to_matrix_row!(
            symbol_ht, ht, UInt32(0), etmp,
            basis.gens[basis.nonred[i]])

        matrix.up2coef[matrix.nrows] = basis.nonred[i]
        # set lead index as 1
        symbol_ht.hashdata[uprows[matrix.nrows][1]].idx = 1
    end

    # needed for correct counting in symbol
    matrix.ncols = matrix.nrows
    matrix.nup = matrix.nrows

    #=
    @warn "after multiplied_poly_to_matrix_row"
    dump(basis, maxdepth=5)
    @warn "matrix"
    dump(matrix, maxdepth=5)
    @warn "ht"
    println(ht.exponents[1:10])
    =#

    symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
    # all pivots are unknown
    for i in symbol_ht.offset:symbol_ht.load
        symbol_ht.hashdata[i].idx = 1
    end

    #=
    @warn "after symbolic"
    dump(basis, maxdepth=5)
    @warn "matrix"
    dump(matrix, maxdepth=5)
    @warn "ht"
    println(ht.exponents[1:10])
    =#

    # x1*x2^6 + 413612941*x1*x2^3
    #  x1^2*x2^3 + 174101409*x1*x2^5
    convert_hashes_to_columns!(matrix, symbol_ht)
    matrix.ncols = matrix.nleft + matrix.nright

    #=
    @warn "after convert"
    dump(matrix, maxdepth=5)
    =#

    sort_matrix_rows_decreasing!(matrix)

    interreduce_matrix_rows!(matrix, basis)

    #=
    @warn "after interreduce"
    dump(matrix, maxdepth=5)
    =#
    # TODO
    # convert_matrix_rows_to_basis_elements_use_symbol!(matrix, basis)
    convert_matrix_rows_to_basis_elements!(matrix, basis, ht, symbol_ht)
    # no longer need in two hashtables
    # TODO
    # ht = symbol_ht

    basis.ntotal = matrix.npivots + basis.ndone
    basis.ndone = matrix.npivots

    #=
    @warn "basis"
    dump(basis, maxdepth=5)
    =#

    #= we may have added some multiples of reduced basis polynomials
    * from the matrix, so we get rid of them. =#
    k = 0
    i = 1
    @label Letsgo
    while i <= basis.ndone
        @inbounds for j in 1:k
            if is_monom_divisible(
                basis.gens[basis.ntotal-i+1][1],
                basis.gens[basis.nonred[j]][1],
                ht)

                i += 1
                @goto Letsgo
            end
        end
        k += 1
        basis.nonred[k] = basis.ntotal - i + 1
        # xd
        basis.lead[k] = ht.hashdata[basis.gens[basis.nonred[k]][1]].divmask
        i += 1
    end
    basis.nlead = k

    # TODO
    # sort_gens_by_lead_increasing_in_reduce!(basis, ht)
end

function select_tobereduced!(
    basis::Basis, tobereduced::Basis,
    matrix::MacaulayMatrix,
    symbol_ht::MonomialHashtable, ht::MonomialHashtable)

    # prepare to load all elems from tobereduced
    # into low rows
    reinitialize_matrix!(matrix, max(basis.ntotal, tobereduced.ntotal))
    resize!(matrix.lowrows, tobereduced.ntotal)

    # TODO
    etmp = zeros(UInt16, ht.explen)

    for i in 1:tobereduced.ntotal
        matrix.nrows += 1
        gen = tobereduced.gens[i]
        h = UInt32(0)
        matrix.lowrows[matrix.nrows] = multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, gen)
        matrix.low2coef[matrix.nrows] = i
    end

    basis.ntotal
    basis.nlead = basis.ndone = basis.ntotal
    basis.isred .= 0
    for i in 1:basis.nlead
        basis.nonred[i] = i
        basis.lead[i] = ht.hashdata[basis.gens[i][1]].divmask
    end
end

#------------------------------------------------------------------------------

function find_multiplied_reducer!(
    basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable,
    vidx::Int)

    e = symbol_ht.exponents[vidx]
    etmp = ht.exponents[1]
    divmask = symbol_ht.hashdata[vidx].divmask

    blen = basis.ndone     #
    leaddiv = basis.lead  #

    # searching for a poly from state.basis whose leading monom
    # divides the given exponent e
    i = 1
    @label Letsgo
    @inbounds while i <= basis.nlead && (leaddiv[i] & ~divmask) != 0
        i += 1
    end

    #@inbounds while i <= basis.nlead && !divv(e, ht.exponents[basis.gens[basis.nonred[i]][1]])
    #    i += 1
    #end

    # here found polynomial from basis with leading monom
    # dividing symbol_ht.exponents[vidx]
    if i <= basis.nlead
        # reducers index and exponent in hash table
        @inbounds rpoly = basis.gens[basis.nonred[i]]
        @inbounds rexp = ht.exponents[rpoly[1]]

        @inbounds for j in 1:ht.explen
            # if it actually does not divide and divmask lies
            if e[j] < rexp[j]
                i += 1
                @goto Letsgo
            end
            etmp[j] = e[j] - rexp[j]
        end
        # now etmp = e // rexp in terms of monomias,
        # hash is linear
        h = symbol_ht.hashdata[vidx].hash - ht.hashdata[rpoly[1]].hash

        matrix.uprows[matrix.nup+1] = multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, rpoly)
        matrix.up2coef[matrix.nup+1] = basis.nonred[i]

        # upsize matrix
        symbol_ht.hashdata[vidx].idx = 2
        matrix.nup += 1
        i += 1
    end

end

#------------------------------------------------------------------------------

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction by G,
# and returns it
function symbolic_preprocessing!(
    basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    # check matrix sizes (we want to omit this I guess)

    symbol_load = symbol_ht.load

    nrr = matrix.ncols
    onrr = matrix.ncols

    #=
    TODO: matrix_enlarge!
    =#
    while matrix.size <= nrr + symbol_load
        matrix.size *= 2
        resize!(matrix.uprows, matrix.size)
        resize!(matrix.up2coef, matrix.size)
    end

    # for each lcm present in symbolic_ht set on select stage
    i = symbol_ht.offset
    #= First round, we add multiplied polynomials which divide  =#
    #= a monomial exponent from selected spairs                 =#
    @inbounds while i <= symbol_load
        # not a reducer already
        if symbol_ht.hashdata[i].idx == 0
            symbol_ht.hashdata[i].idx = 1
            matrix.ncols += 1
            find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        end
        i += 1
    end

    #= Second round, we add multiplied polynomials which divide  =#
    #= lcm added on previous for loop                            =#
    while i <= symbol_ht.load
        if matrix.size == matrix.nup
            matrix.size *= 2
            # TODO:
            resize!(matrix.uprows, matrix.size)
            resize!(matrix.up2coef, matrix.size)
        end

        symbol_ht.hashdata[i].idx = 1
        matrix.ncols += 1
        find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        i += 1
    end

    # shrink matrix sizes, set constants
    resize!(matrix.uprows, matrix.nup)
    # resize!(matrix.up2coef, matrix.nup)

    matrix.nrows += matrix.nup - onrr
    matrix.nlow = matrix.nrows - matrix.nup
    matrix.size = matrix.nrows

end

#------------------------------------------------------------------------------

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction by G,
# and returns it
function symbolic_preprocessing_relaxed!(
    basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    # check matrix sizes (we want to omit this I guess)

    symbol_load = symbol_ht.load

    nrr = matrix.ncols
    onrr = matrix.ncols

    #=
    TODO: matrix_enlarge!
    =#
    while matrix.size <= nrr + symbol_load
        matrix.size *= 2
        resize!(matrix.uprows, matrix.size)
        resize!(matrix.up2coef, matrix.size)
    end

    # for each lcm present in symbolic_ht set on select stage
    i = symbol_ht.offset
    #= First round, we add multiplied polynomials which divide  =#
    #= a monomial exponent from selected spairs                 =#
    @inbounds while i <= symbol_load
        # not a reducer already
        if symbol_ht.hashdata[i].idx == 0
            symbol_ht.hashdata[i].idx = 1
            matrix.ncols += 1
            find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
        end
        i += 1
    end

    #= Second round, we add multiplied polynomials which divide  =#
    #= lcm added on previous for loop                             =#
    #=
    while i <= symbol_ht.load
    if matrix.size == matrix.nup
        matrix.size *= 2
        resize!(matrix.uprows, matrix.size)
        resize!(matrix.up2coef, matrix.size)
    end

    symbol_ht.hashdata[i].idx = 1
    matrix.ncols += 1
    find_multiplied_reducer!(basis, matrix, ht, symbol_ht, i)
    i += 1
    end
    =#

    # shrink matrix sizes, set constants
    resize!(matrix.uprows, matrix.nup)
    # resize!(matrix.up2coef, matrix.nup)

    matrix.nrows += matrix.nup - onrr
    matrix.nlow = matrix.nrows - matrix.nup
    matrix.size = matrix.nrows

end


function select_normal_and_discard!(pairset::Pairset, basis::Basis,
    matrix::MacaulayMatrix, ht::MonomialHashtable,
    symbol_ht::MonomialHashtable; maxpairs::Int=0)

    sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
    # sort by degree

    ps = pairset.pairs
    min_deg = ps[1].deg
    min_idx = 0

    while min_idx < pairset.load && ps[min_idx+1].deg == min_deg
        min_idx += 1
    end

    # number of discarded pairs
    npairs = min_idx

    @debug "Discarder $(npairs) pairs"

    for i in 1:pairset.load-npairs
        ps[i] = ps[i+npairs]
    end
    pairset.load -= npairs
end

function select_normal!(
    pairset::Pairset, basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable;
    maxpairs::Int=0)

    sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
    # sort by degree

    ps = pairset.pairs
    min_deg = ps[1].deg
    min_idx = 0

    while min_idx < pairset.load && ps[min_idx+1].deg == min_deg
        min_idx += 1
    end

    # number of selected pairs
    # npairs = min(min_idx, 20)
    npairs = min_idx

    @debug "Selected $(npairs) pairs"

    sort_pairset_by_lcm!(pairset, npairs, ht)

    reinitialize_matrix!(matrix, npairs)

    uprows = matrix.uprows
    lowrows = matrix.lowrows

    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens = Vector{Int}(undef, 2 * npairs)

    # buffer !
    etmp = ht.exponents[1]
    i = 1
    while i <= npairs
        matrix.ncols += 1
        load = 1
        lcm = ps[i].lcm
        j = i

        # we collect all generators with same lcm into gens
        @inbounds while j <= npairs && ps[j].lcm == lcm
            gens[load] = ps[j].poly1
            load += 1
            gens[load] = ps[j].poly2
            load += 1
            j += 1
        end
        load -= 1

        # sort by number in the basis (by=identity)
        sort_generators_by_position!(gens, load)

        # now we collect reducers, and reduced

        # first generator index in groebner basis
        prev = gens[1]
        # first generator in hash table
        poly = basis.gens[prev]
        # first generator lead monomial index in hash data
        vidx = poly[1]

        # first generator exponent
        eidx = ht.exponents[vidx]
        # exponent of lcm corresponding to first generator
        elcm = ht.exponents[lcm]
        for u in 1:ht.explen
            etmp[u] = elcm[u] - eidx[u]
        end
        # now etmp contents complement to eidx in elcm

        # hash of complement
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # add row as a reducer
        matrix.nup += 1
        uprows[matrix.nup] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
        # map upper row to index in basis
        matrix.up2coef[matrix.nup] = prev

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
            elcm = ht.exponents[lcm]

            # index in gb
            prev = gens[k]
            # poly of indices of monoms in hash table
            poly = basis.gens[prev]
            vidx = poly[1]
            # leading monom idx
            eidx = ht.exponents[vidx]
            for u in 1:ht.explen
                etmp[u] = elcm[u] - eidx[u]
            end

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlow += 1
            lowrows[matrix.nlow] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
            # map lower row to index in basis
            matrix.low2coef[matrix.nlow] = prev

            symbol_ht.hashdata[lowrows[matrix.nlow][1]].idx = 2

            matrix.nrows += 1
        end

        i = j
    end

    resize!(matrix.lowrows, matrix.nrows - matrix.ncols)

    # remove selected parirs from pairset
    for i in 1:pairset.load-npairs
        ps[i] = ps[i+npairs]
    end
    pairset.load -= npairs
end

#------------------------------------------------------------------------------

function select_isgroebner!(
    pairset::Pairset, basis::Basis, matrix::MacaulayMatrix,
    ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    ps = pairset.pairs
    npairs = pairset.load

    sort_pairset_by_lcm!(pairset, npairs, ht)

    reinitialize_matrix!(matrix, npairs)

    uprows = matrix.uprows
    lowrows = matrix.lowrows

    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens = Vector{Int}(undef, 2 * npairs)

    # buffer !
    etmp = ht.exponents[1]
    i = 1
    while i <= npairs
        matrix.ncols += 1
        load = 1
        lcm = ps[i].lcm
        j = i

        # we collect all generators with same lcm into gens
        @inbounds while j <= npairs && ps[j].lcm == lcm
            gens[load] = ps[j].poly1
            load += 1
            gens[load] = ps[j].poly2
            load += 1
            j += 1
        end
        load -= 1

        # sort by number in the basis (by=identity)
        sort_generators_by_position!(gens, load)

        # now we collect reducers, and reduced

        # first generator index in groebner basis
        prev = gens[1]
        # first generator in hash table
        poly = basis.gens[prev]
        # first generator lead monomial index in hash data
        vidx = poly[1]

        # first generator exponent
        eidx = ht.exponents[vidx]
        # exponent of lcm corresponding to first generator
        elcm = ht.exponents[lcm]
        for u in 1:ht.explen
            etmp[u] = elcm[u] - eidx[u]
        end
        # now etmp contents complement to eidx in elcm

        # hash of complement
        htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

        # add row as a reducer
        matrix.nup += 1
        uprows[matrix.nup] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
        # map upper row to index in basis
        matrix.up2coef[matrix.nup] = prev

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
            elcm = ht.exponents[lcm]

            # index in gb
            prev = gens[k]
            # poly of indices of monoms in hash table
            poly = basis.gens[prev]
            vidx = poly[1]
            # leading monom idx
            eidx = ht.exponents[vidx]
            for u in 1:ht.explen
                etmp[u] = elcm[u] - eidx[u]
            end

            htmp = ht.hashdata[lcm].hash - ht.hashdata[vidx].hash

            # add row to be reduced
            matrix.nlow += 1
            lowrows[matrix.nlow] = multiplied_poly_to_matrix_row!(symbol_ht, ht, htmp, etmp, poly)
            # map lower row to index in basis
            matrix.low2coef[matrix.nlow] = prev

            symbol_ht.hashdata[lowrows[matrix.nlow][1]].idx = 2

            matrix.nrows += 1
        end

        i = j
    end

    resize!(matrix.lowrows, matrix.nrows - matrix.ncols)
end

#------------------------------------------------------------------------------

#=
    Input ivariants:
        - ring is set, and ring.ch == basis.ch, and ring ~ ht
        - divmasks in ht are set
        - basis is filled so that
            basis.ntotal = actual number of elements
            basis.ndone  = 0
            basis.nlead  = 0

    Output invariants:
        - basis.ndone == basis.ntotal == basis.nlead
        - basis.gens and basis.coeffs of size basis.ndone
        - basis elements are sorted increasingly wrt ordering on lead elements
        - divmasks in basis are filled and coincide to divmasks in hashtable

=#
function f4!(ring::PolyRing,
    basis::Basis{Coefftype},
    tracer::Tracer,
    pairset,
    ht,
    reduced,
    linalg,
    rng) where {Coefftype<:Coeff}

   # print("input: $(basis.ntotal) gens, $(ring.nvars) vars. ")

   @assert ring.ch == basis.ch
   @assert ring.ord == ht.ord && ring.nvars == ht.nvars && ring.explen == ht.explen
   # @error "hashtable divmasks"
   # println(ht.exponents[2:ht.load])
   # println("###########")
   # println(map(x->bitstring(x.divmask), ht.hashdata[2:ht.load]))
   # println("###########")
   @assert basis.ndone == 0

   # matrix storing coefficients in rows
   # wrt columns representing the current monomial basis
   matrix = initialize_matrix(ring, Coefftype)

   # initialize hash tables for update and symbolic preprocessing steps
   update_ht  = initialize_secondary_hash_table(ht)
   symbol_ht  = initialize_secondary_hash_table(ht)

   # a set to store critical pairs of polynomials to be reduced
   # if tracer.ready
   #     pairset_initial_size = tracer.pairset_size
   # else
   #     pairset_initial_size = 64
   # end
   # pairset = initialize_pairset(initial_size=pairset_initial_size)

   # @warn "ht initialily" ht.load ht.size

   # makes basis fields valid,
   # does not copy,
   # checks for redundancy of new elems
   plcm = Vector{ExponentIdx}(undef, 0)
   if tracer.ready
       resize!(plcm, tracer.basis_ntotal + 1)
   end

   # @error "" length(pairset.pairs) length(plcm)

   pairset_size = update!(pairset, basis, ht, update_ht, plcm)
   if !tracer.ready
       tracer.pairset_size = pairset_size
   end

   # @warn "tracer info" tracer.ready tracer.pairset_size

   # @warn "ht update" ht.load ht.size

   d = 0
   # while there are pairs to be reduced
   while !isempty(pairset)
       d += 1
       @debug "F4 ITER $d"
       @debug "Available $(pairset.load) pairs"

       # TODO: learn, and select\discard S-polynomials
       if tracer.ready && tracer.isredundant_iter[d] == 0
           select_normal_and_discard!(pairset, basis, matrix, ht, symbol_ht)
           matrix    = initialize_matrix(ring, Coefftype)
           symbol_ht = initialize_secondary_hash_table(ht)
           continue
       end

       # selects pairs for reduction from pairset following normal strategy
       # (minimal lcm degrees are selected),
       # and puts these into the matrix rows
       select_normal!(pairset, basis, matrix, ht, symbol_ht)

       # @warn "ht select" ht.load ht.size

       symbolic_preprocessing!(basis, matrix, ht, symbol_ht)
       # symbolic_preprocessing_relaxed!(basis, matrix, ht, symbol_ht)
       @debug "Matrix of size $((matrix.nrows, matrix.ncols)), density TODO"

       # @warn "ht symbolic" ht.load ht.size

       # reduces polys and obtains new potential basis elements
       reduction!(basis, matrix, ht, symbol_ht, linalg, rng)
       @debug "Matrix reduced, density TODO"

       if !tracer.ready
           push!(tracer.isredundant_iter, 0)
           if matrix.npivots != 0
               tracer.isredundant_iter[end] = 1
           end
       end

       #=
       if tracer.ready && length(tracer.isredundant_iter) > d
           if tracer.isredundant_iter[d + 1] == 0
               matrix    = initialize_matrix(ring, Coefftype)
               symbol_ht = initialize_secondary_hash_table(ht)
               continue
           end
       end
       =#

       # update the current basis with polynomials produced from reduction,
       # does not copy,
       # checks for redundancy
       pairset_size = update!(pairset, basis, ht, update_ht, plcm)
       if !tracer.ready
           tracer.pairset_size = max(pairset_size, tracer.pairset_size)
       end
       # @warn "ht update" ht.load ht.size

       # TODO: is this okay hm ?
       # to be changed
       # TODO: clean hashtable
       matrix    = initialize_matrix(ring, Coefftype)
       symbol_ht = initialize_secondary_hash_table(ht)
       # clear symbolic hashtable
       # clear matrix

       if d > 10000
           @error "Something is probably wrong in f4. Please submit an issue."
           break
       end
   end

   tracer.ready = true
   tracer.basis_ntotal = basis.ntotal

   # remove redundant elements
   filter_redundant!(basis)

   if reduced
       reducegb_f4!(basis, matrix, ht, symbol_ht)
   end

   standardize_basis!(basis, ht, ht.ord)

   # assertion
   #=
   println("#####")
   println(basis.lead)
   println(map(x->bitstring(x.divmask), ht.hashdata[2:ht.load]))
   println("#####")
   =#
   #=
   for i in 1:basis.nlead
       #=
       println(basis.lead[i])
       println(basis.gens[i][1])
       println(ht.exponents[basis.gens[i][1]])
       println(UInt32(ht.hashdata[basis.gens[i][1]].divmask))
       =#
       @assert basis.lead[i] == ht.hashdata[basis.gens[i][1]].divmask
   end
   =#
   # println("HT: $(ht.load)/$(ht.size)")

   nothing
end

function f4!(ring::PolyRing,
        basis::Basis{Coefftype},
        ht,
        reduced,
        linalg,
        rng) where {Coefftype}
    ps = initialize_pairset()
    tracer = Tracer()
    f4!(ring, basis, tracer, ps, ht, reduced, linalg, rng)
end
