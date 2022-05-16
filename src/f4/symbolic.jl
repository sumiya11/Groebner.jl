
function select_normal!(
            pairset::Pairset, basis::Basis, matrix::MacaulayMatrix,
            ht::MonomialHashtable, symbol_ht::MonomialHashtable;
            maxpairs::Int=0)

    sort_pairset_by_degree!(pairset, 1, pairset.load - 1)
    # sort by degree

    ps = pairset.pairs
    min_deg = ps[1].deg
    min_idx = 0

    while min_idx < pairset.load && ps[min_idx + 1].deg == min_deg
        min_idx += 1
    end

    # number of selected pairs
    npairs = min_idx

    sort_pairset_by_lcm!(pairset, npairs, ht)

    reinitialize_matrix!(matrix, npairs)

    uprows  = matrix.uprows
    lowrows = matrix.lowrows

    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens   = Vector{Int}(undef, 2*npairs)

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
    for i in 1:pairset.load - npairs
        ps[i] = ps[i + npairs]
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

    uprows  = matrix.uprows
    lowrows = matrix.lowrows

    # polynomials from pairs in order (p11, p12)(p21, p21)
    # (future rows of the matrix)
    gens   = Vector{Int}(undef, 2*npairs)

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
        @inbounds rexp  = ht.exponents[rpoly[1]]

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

        matrix.uprows[matrix.nup + 1] = multiplied_poly_to_matrix_row!(symbol_ht, ht, h, etmp, rpoly)
        matrix.up2coef[matrix.nup + 1] = basis.nonred[i]

        # upsize matrix
        symbol_ht.hashdata[vidx].idx = 2
        matrix.nup += 1
        i += 1
    end

end

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction by G,
# and returns it
function symbolic_preprocessing!(
                basis::Basis, matrix::MacaulayMatrix,
                ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    # check matrix sizes (we want to omit this I guess)

    symbol_load = symbol_ht.load

    nrr    = matrix.ncols
    onrr   = matrix.ncols

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

    # shrink matrix sizes, set constants
    resize!(matrix.uprows, matrix.nup)
    # resize!(matrix.up2coef, matrix.nup)

    matrix.nrows += matrix.nup - onrr
    matrix.nlow  = matrix.nrows - matrix.nup
    matrix.size  = matrix.nrows

end

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction by G,
# and returns it
function symbolic_preprocessing_relaxed!(
                basis::Basis, matrix::MacaulayMatrix,
                ht::MonomialHashtable, symbol_ht::MonomialHashtable)

    # check matrix sizes (we want to omit this I guess)

    symbol_load = symbol_ht.load

    nrr    = matrix.ncols
    onrr   = matrix.ncols

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
    matrix.nlow  = matrix.nrows - matrix.nup
    matrix.size  = matrix.nrows

end

function update_pairset!(
            pairset,
            basis,
            ht,
            update_ht,
            idx)

    pl = pairset.load
    bl = idx
    nl = pl + bl
    ps = pairset.pairs

    new_lead  = basis.gens[idx][1]

    # initialize new critical lcms
    plcm = Vector{Int}(undef, bl + 1)

    # for each combination (new_Lead, basis.gens[i][1])
    # generate a pair
    @inbounds for i in 1:bl-1
        plcm[i] = get_lcm(basis.gens[i][1], new_lead, ht, update_ht)
        deg = update_ht.hashdata[plcm[i]].deg
        newidx = pl + i
        if basis.isred[i] == 0
            ps[newidx] = SPair(i, idx, plcm[i], deg)
        else
            # lcm == 0 will mark redundancy of spair
            ps[newidx] = SPair(i, idx, 0, deg)
        end
    end

    # traverse existing pairs
    @inbounds for i in 1:pl
        j = ps[i].poly1
        l = ps[i].poly2
        m = max(ps[pl + l].deg, ps[pl + j].deg)

        # if an existing pair is divisible by lead of new poly
        # and has a greater degree than newly generated one
        if is_monom_divisible(ps[i].lcm, new_lead, ht) && ps[i].deg > m
            # mark lcm as 0
            ps[i] = SPair(ps[i].poly1, ps[i].poly2, 0, ps[i].deg)
        end
    end
    # TODO: this can be done faster is we
    # create only non-redundant pairs at first

    # traverse new pairs to check for redundancy
    j = 1
    for i in 1:bl-1
        if basis.isred[i] == 0
            ps[pl + j] = ps[pl + i]
            j += 1
        end
    end

    sort_pairset_by_degree!(pairset, pl + 1, j - 2)

    @inbounds for i in 1:j - 1
        plcm[i] = ps[pl + i].lcm
    end
    plcm[j] = 0
    pc = j
    pc -= 1

    # mark redundancy of some pairs from plcm
    @inbounds for j in 1:pc
        # if is not redundant
        if plcm[j] > 0
            check_monomial_division_in_update(plcm, j + 1, pc, plcm[j], update_ht)
        end
    end

    # remove useless pairs from pairset
    # by moving them to the end
    j = 1
    @inbounds for i in 1:pairset.load
        ps[i].lcm == 0 && continue
        ps[j] = ps[i]
        j += 1
    end

    # assure that basis hashtable can store new lcms
    if ht.size - ht.load <= pc
        enlarge_hash_table!(ht)
    end

    # add new lcms to the basis hashtable,
    # including j and not including pc
    insert_plcms_in_basis_hash_table!(pairset, pl, ht, update_ht, basis, plcm, j, pc+1)

    # mark redundant elements in masis
    nonred = basis.nonred
    lml = basis.nlead
    @inbounds for i in 1:lml
        if basis.isred[nonred[i]] == 0
            if is_monom_divisible(basis.gens[nonred[i]][1], new_lead, ht)
                basis.isred[nonred[i]] = 1
            end
        end
    end

end

function update_basis!(
            basis,
            ht::MonomialHashtable,
            update_ht::MonomialHashtable)

    # here we could check overall redundancy and update basis.lead

    k = 1
    lead   = basis.lead
    nonred = basis.nonred

    for i in 1:basis.nlead
        if basis.isred[nonred[i]] == 0
            basis.lead[k]   = lead[i]
            basis.nonred[k] = nonred[i]
            k += 1
        end
    end
    basis.nlead = k - 1

    for i in basis.ndone+1:basis.ntotal
        if basis.isred[i] == 0
            lead[k]   = ht.hashdata[basis.gens[i][1]].divmask
            nonred[k] = i
            k += 1
        end
    end

    basis.nlead = k - 1
    basis.ndone = basis.ntotal
end

# checks if element of basis at position idx is redundant
function is_redundant!(
            pairset, basis, ht, update_ht, idx)


    reinitialize_hash_table!(update_ht, 2*idx)

    ps = pairset.pairs

    # lead of new polynomial
    lead_new = basis.gens[idx][1]
    # degree of lead
    lead_deg = ht.hashdata[lead_new].deg

    for i in idx+1:basis.ntotal
        i == idx && continue
        if basis.isred[i] == 1
            continue
        end

        lead_i = basis.gens[i][1]

        if is_monom_divisible(lead_new, lead_i, ht)
            lcm_new = get_lcm(lead_i, lead_new, ht, ht)

            psidx = pairset.load + 1
            ps[psidx] = SPair(i, idx, lcm_new, ht.hashdata[lcm_new].deg)

            basis.isred[idx] = 1
            pairset.load += 1

            return true
        end
    end

    return false
end

function update!(
        pairset::Pairset,
        basis::Basis,
        ht::MonomialHashtable,
        update_ht::MonomialHashtable)

    #=
        Always check redundancy, for now
    =#

    # number of added elements
    npivs = basis.ntotal

    # number of potential critical pairs to add
    npairs = basis.ndone * npivs
    for i in 1:npivs
        npairs += i
    end


    # make sure pairset and update hashtable have enough
    # space to store new pairs
    # TODO: we create too big array, can be fixed
    check_enlarge_pairset!(pairset, npairs)

    # for each new element in basis
    for i in basis.ndone+1:basis.ntotal
        # check redundancy of new poly
        if is_redundant!(pairset, basis, ht, update_ht, i)
            continue
        end
        update_pairset!(pairset, basis, ht, update_ht, i)
    end

    update_basis!(basis, ht, update_ht)
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
