
#------------------------------------------------------------------------------

function insert_in_hash_table!(ht::MonomialHashtable, e::ExponentVector)
    # generate hash
    he = UInt32(0)

    # here, e[i] is of type UInt16, while hasher[i] is UInt32 =(
    @inbounds for i in 1:ht.explen
        he += ht.hasher[i] * e[i]
    end

    # find new elem position in the table
    hidx = Int(he)  # Int for type stability
    # power of twoooo
    mod = UInt32(ht.size - 1)
    i = UInt32(1)

    @label Restart
    while i < ht.size
        hidx = hashnextindex(he, i, mod)
        @inbounds vidx = ht.hashtable[hidx]

        # if free
        vidx == 0 && break

        # if not free and not same hash
        if ht.hashdata[vidx].hash != he
            i += UInt32(1)
            continue
        end

        present = ht.exponents[vidx]
        @inbounds for j in 1:ht.explen
            # if hash collision
            if present[j] != e[j]
                i += UInt32(1)
                @goto Restart
            end
        end

        # already present in hashtable
        return vidx
    end

    # add its position to hashtable, and insert exponent to that position
    vidx = ht.load + 1
    ht.hashtable[hidx] = vidx

    # TODO: check efficiency
    ht.exponents[vidx] = similar(e)
    ve = ht.exponents[vidx]
    @inbounds for i in 1:length(e)
        ve[i] = e[i]
    end
    # ht.exponents[vidx] = copy(e)

    divmask = generate_monomial_divmask(e, ht)
    ht.hashdata[vidx] = Hashvalue(he, divmask, 0, e[end])

    ht.load += 1

    return vidx
end

#------------------------------------------------------------------------------

#=
    Having `ht` filled with monomials from input polys,
    computes ht.divmap and divmask for each of already stored monomials
=#
function fill_divmask!(ht::MonomialHashtable)
    ndivvars = ht.ndivvars
    divvars = ht.divvars

    min_exp = Vector{UInt16}(undef, ndivvars)
    max_exp = Vector{UInt16}(undef, ndivvars)

    e = ht.exponents[ht.offset]
    for i in 1:ndivvars
        min_exp[i] = e[divvars[i]]
        max_exp[i] = e[divvars[i]]
    end

    for i in ht.offset:ht.load # TODO: offset
        e = ht.exponents[i]
        for j in 1:ndivvars
            if e[divvars[j]] > max_exp[j]
                max_exp[j] = e[divvars[j]]
                continue
            end
            if e[divvars[j]] < min_exp[j]
                min_exp[j] = e[divvars[j]]
            end
        end
    end

    ctr = 1
    steps = UInt32(0)
    for i in 1:ndivvars
        steps = div(max_exp[i] - min_exp[i], ht.ndivbits)
        (steps == 0) && (steps += 1)
        for j in 1:ht.ndivbits
            ht.divmap[ctr] = steps
            steps += 1
            ctr += 1
        end
    end

    for vidx in ht.offset:ht.load
        unmasked = ht.hashdata[vidx]
        e = ht.exponents[vidx]
        divmask = generate_monomial_divmask(e, ht)
        ht.hashdata[vidx] = Hashvalue(unmasked.hash, divmask, 0, e[end])
    end
end

#=
    TODO

=#
function generate_monomial_divmask(
    e::Vector{UInt16},
    ht::MonomialHashtable)

    divvars = ht.divvars
    divmap = ht.divmap

    ctr = UInt32(1)
    res = UInt32(0)
    for i in 1:ht.ndivvars
        for j in 1:ht.ndivbits
            @inbounds if e[divvars[i]] >= divmap[ctr]
                res |= UInt32(1) << (ctr - 1)
            end
            ctr += UInt32(1)  # for type stability
        end
    end

    res
end

#------------------------------------------------------------------------------

# h1 divisible by h2
function is_monom_divisible(h1::Int, h2::Int, ht::MonomialHashtable)

    @inbounds if (ht.hashdata[h2].divmask & ~ht.hashdata[h1].divmask) != 0
        return false
    end

    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]
    # TODO: one less iteration is possible
    @inbounds for i in 1:ht.explen
        if e1[i] < e2[i]
            return false
        end
    end

    return true
end

function is_gcd_const(h1::Int, h2::Int, ht::MonomialHashtable)
    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]

    for i in 1:ht.explen-1
        if e1[i] != 0 && e2[i] != 0
            return false
        end
    end

    return true
end

#------------------------------------------------------------------------------

# compare pairwise divisibility of lcms from a[first:last] with lcm
function check_monomial_division_in_update(
    a::Vector{Int}, first::Int, last::Int,
    lcm::Int, ht::MonomialHashtable)

    # pairs are sorted, we only need to check entries above starting point

    divmask = ht.hashdata[lcm].divmask
    lcmexp = ht.exponents[lcm]

    j = first
    @label Restart
    while j <= last
        # bad lcm
        if a[j] == 0
            j += 1
            continue
        end
        # fast division check
        @inbounds if (~ht.hashdata[a[j]].divmask & divmask) != 0
            j += 1
            continue
        end
        @inbounds ea = ht.exponents[a[j]]
        @inbounds for i in 1:ht.explen
            if ea[i] < lcmexp[i]
                j += 1
                @goto Restart
            end
        end
        # mark as redundant
        a[j] = 0
    end

end

#------------------------------------------------------------------------------

# add monomials from `poly` multiplied by exponent vector `etmp`
# with hash `htmp` to hashtable `symbol_ht`,
# and substitute hashes in row
function insert_multiplied_poly_in_hash_table!(
    row::Vector{Int},
    htmp::UInt32,
    etmp::ExponentVector,
    poly::Vector{Int},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable)

    # oof

    # length of poly to add
    len = length(poly)
    explen = ht.explen

    mod = UInt32(symbol_ht.size - 1)

    bexps = ht.exponents
    bdata = ht.hashdata

    sexps = symbol_ht.exponents
    sdata = symbol_ht.hashdata

    # @error "" symbol_ht.load symbol_ht.size symbol_ht.load / symbol_ht.size

    l = 1 # hardcoding 1 does not seem nice =(
    @label Letsgo
    while l <= len
        # we iterate over all monoms of the given poly,
        # multiplying them by htmp/etmp,
        # and inserting into symbolic hashtable

        # hash is linear, so that
        # hash(e1 + e2) = hash(e1) + hash(e2)
        # We also assume that the hashing vector is shared same
        # between all created hashtables
        @inbounds h = htmp + bdata[poly[l]].hash
        # TODO! -- check mult hash in the table

        @inbounds e = bexps[poly[l]]
        # println("monom of index $(poly[l]) : $e")

        lastidx = symbol_ht.load + 1
        #=
        if !isassigned(sexps, lastidx)
            sexps[lastidx] = Vector{UInt16}(undef, explen)
        end
        enew = sexps[lastidx]
        =#
        @inbounds enew = sexps[1]

        @inbounds for j in 1:explen
            # multiplied monom exponent
            enew[j] = etmp[j] + e[j]
        end

        # now insert into hashtable
        k = h

        i = UInt32(1)

        @label Restart
        while i <= symbol_ht.size  # TODO: < or <= ?
            k = hashnextindex(h, i, mod)

            @inbounds vidx = symbol_ht.hashtable[k]
            # if index is free
            vidx == 0 && break
            # if different exponent is stored here
            @inbounds if sdata[vidx].hash != h

                # global ADD_ROW_COLLISION
                # ADD_ROW_COLLISION += 1

                i += UInt32(1)
                continue
            end

            @inbounds estored = sexps[vidx]
            @inbounds for j in 1:explen
                # hash collision, restarting search
                if estored[j] != enew[j]
                    i += UInt32(1)

                    # global ADD_ROW_COLLISION
                    # ADD_ROW_COLLISION += 1

                    @goto Restart
                end
            end

            # @error "hit"

            # global ADD_ROW_HIT
            # ADD_ROW_HIT += 1

            @inbounds row[l] = vidx
            l += 1

            @goto Letsgo
        end

        # global ADD_ROW_MISS
        # ADD_ROW_MISS += 1

        # @warn "miss"

        # add multiplied exponent to hash table
        if !isassigned(sexps, lastidx)
            sexps[lastidx] = ExponentVector(undef, explen)
        end
        sexpsnew = sexps[lastidx]
        for j in 1:explen
            # multiplied monom exponent
            @inbounds sexpsnew[j] = enew[j]
        end
        symbol_ht.hashtable[k] = lastidx

        divmask = generate_monomial_divmask(enew, symbol_ht)
        sdata[lastidx] = Hashvalue(h, divmask, 0, enew[end])

        row[l] = lastidx
        l += 1
        symbol_ht.load += 1
    end

    row
end

# TODO: exponent vectors and coeffs everywhere

function multiplied_poly_to_matrix_row!(
    symbolic_ht::MonomialHashtable, basis_ht::MonomialHashtable,
    htmp::UInt32, etmp::Vector{UInt16}, poly::Vector{Int})

    row = similar(poly)
    while 1.4 * (symbolic_ht.load + length(poly)) >= symbolic_ht.size
        enlarge_hash_table!(symbolic_ht)
    end

    #=
    load = symbolic_ht.load
    size = symbolic_ht.size
    if load % 100 == 0 && load/size > 0.7
        @warn "inserting in a very dense table" load size load/size
    end
    =#

    insert_multiplied_poly_in_hash_table!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
    #if load % 100 == 0
    #    println("time: ", time() - t1)
    #end
end

#------------------------------------------------------------------------------

function insert_in_basis_hash_table_pivots(
    row::Vector{Int},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    col2hash::Vector{Int})

    while ht.size - ht.load <= length(row)
        enlarge_hash_table!(ht)
    end

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.exponents

    mod = UInt32(ht.size - 1)
    explen = ht.explen
    bdata = ht.hashdata
    bexps = ht.exponents
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    while l <= length(row)
        hidx = col2hash[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        # TODO: speed this up
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = UInt32(1)
        @label Restart
        while i <= ht.size
            k = hashnextindex(h, i, mod)
            hm = bhash[k]

            hm == 0 && break
            @inbounds if bdata[hm].hash != h
                i += UInt32(1)
                continue
            end

            ehm = bexps[hm]
            @inbounds for j in 1:explen
                if e[j] != ehm[j]
                    i += UInt32(1)
                    @goto Restart
                end
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(h, sdata[hidx].divmask,
            sdata[hidx].idx, sdata[hidx].deg)

        ht.load += 1
    end
end

function insert_plcms_in_basis_hash_table!(
    pairset::Pairset,
    off::Int,
    ht::MonomialHashtable,
    update_ht::MonomialHashtable,
    basis::Basis,
    plcm::Vector{Int},
    ifirst::Int, ilast::Int)

    # including ifirst and not including ilast

    gens = basis.gens
    mod = UInt32(ht.size - 1)
    ps = pairset.pairs

    m = ifirst
    l = 1
    @label Letsgo
    while l < ilast
        if plcm[l] == 0
            l += 1
            continue
        end

        if is_gcd_const(gens[ps[off+l].poly1][1], gens[ps[off+1].poly2][1], ht)
            l += 1
            continue
        end

        ps[m] = ps[off+l]

        # TODO: IT IS NOT CORRECT
        # upd: it is, but it can be done better
        h = update_ht.hashdata[plcm[l]].hash
        ht.exponents[ht.load+1] = copy(update_ht.exponents[plcm[l]])
        n = ht.exponents[ht.load+1]

        k = h
        i = UInt32(1)
        @label Restart
        while i <= ht.size
            k = hashnextindex(h, i, mod)
            hm = ht.hashtable[k]

            hm == 0 && break
            if ht.hashdata[hm].hash != h
                i += UInt32(1)
                continue
            end

            ehm = ht.exponents[hm]

            # @info "SO, we have " n ehm m
            for j in 1:ht.explen
                if ehm[j] != n[j]
                    i += UInt32(1)
                    @goto Restart
                end
            end

            ps[m] = SPair(ps[m].poly1, ps[m].poly2, hm, ps[m].deg)
            m += 1
            l += 1
            @goto Letsgo
        end

        ht.hashtable[k] = pos = ht.load + 1

        uhd = update_ht.hashdata

        ll = plcm[l]
        ht.hashdata[ht.load+1] = Hashvalue(h, uhd[ll].divmask, 0, uhd[ll].deg)

        ht.load += 1
        ps[m] = SPair(ps[m].poly1, ps[m].poly2, pos, ps[m].deg)
        m += 1
        l += 1
    end

    pairset.load = m - 1
end

#------------------------------------------------------------------------------

# computes lcm of he1 and he2 as exponent vectors from ht1
# and inserts it in ht2
function get_lcm(he1::Int, he2::Int,
    ht1::MonomialHashtable, ht2::MonomialHashtable)

    @inbounds e1 = ht1.exponents[he1]
    @inbounds e2 = ht1.exponents[he2]
    @inbounds etmp = ht1.exponents[1]

    # TODO: degrevlex only
    @inbounds etmp[end] = 0
    @inbounds for i in 1:ht1.explen-1
        etmp[i] = max(e1[i], e2[i])
        etmp[end] += etmp[i]
    end

    insert_in_hash_table!(ht2, etmp)
end

#------------------------------------------------------------------------------
