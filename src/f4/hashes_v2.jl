

# TODO: make immutable
#=
    stores hash of a monomial,
    corresponding divmask to speed up divisibility checks,
    index for position matrix (defaults to zero),
    and the todal degree
=#
mutable struct Hashvalue
    hash::UInt32

    #=

    =#
    divmask::UInt32

    idx::Int
    deg::Int
end

mutable struct CustomMonomialHashtable
    # stores monomial exponent vectors,
    # assumes degrees are < 2^16
    exponents::Vector{Vector{UInt16}}

    # maps exponent hash to its position in exponents array
    hashtable::Vector{Int}

    # stores hashes, division masks,
    # and other valuable info
    # for each hashtable enrty
    hashdata::Vector{Hashvalue}

    # values to hash exponents with, i.e
    # hash(e) = sum(hasher .* e)
    hasher::Vector{UInt32}

    #= Ring information =#
    # number of variables
    nvars::Int
    # raw length of exponent vector
    explen::Int
    # ring monomial ordering
    ord::Symbol

    #= Divisibility =#
    # divisor map to check divisibility faster
    divmap::Vector{UInt32}
    # variables selected for divmap
    divvars::Vector{Int}
    # count of divmap variables
    ndivvars::Int
    # bits per div variable
    ndivbits::Int

    size::Int
    # elements added
    load::Int
    #
    offset::Int
end

#------------------------------------------------------------------------------

# initialize and set fields for basis hashtable
# TODO: Initial size is a nice parameter to play with
# TODO: Get rid of `Tv` dependency
function initialize_basis_hash_table(ring::PolyRing{Tv}; initial_size=2^3) where {Tv}
    exponents = Vector{Vector{UInt16}}(undef, initial_size)
    hashdata  = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(Int, initial_size)

    nvars = ring.nvars
    explen = ring.explen
    ord = ring.ord

    # initialize hashing vector
    hasher = zeros(UInt32, explen)
    for i in 1:explen
        # we don't want hash vector components to be zero
        while hasher[i] == 0
            hasher[i] = rand(UInt32)
        end
    end

    # exponents[1:load] cover all stored exponents
    # , also exponents[1] is zeroed by default
    load = 1
    size = initial_size

    # exponents array starts from index offset,
    # We store buffer array at index 1
    offset = 2

    # initialize fast divisibility params
    charbit   = 8 # TODO ??
    int32bits = charbit * sizeof(Int32)
    int32bits!= 32 && error("Strange story with ints")
    ndivbits  = div(int32bits, nvars)
    # division mask stores at least 1 bit
    # per each of first charbit*sizeof(Int32) variables
    ndivbits == 0 && (ndivbits += 1)
    ndivvars = nvars < int32bits ? nvars : int32bits
    divvars  = Vector{Int}(undef, ndivvars)
    divmap   = Vector{UInt32}(undef, ndivvars * ndivbits)
    # count only first ndivvars variables for divisibility checks
    for i in 1:ndivvars
        divvars[i] = i
    end

    # first stored exponent used as buffer lately
    exponents[1] = zeros(UInt16, explen)

    CustomMonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

# initialize hashtable either for `symbolic_preprocessing` or for `update` functions
# These are of the same purpose and structure as basis hashtable,
# but are more local oriented
function initialize_secondary_hash_table(basis_ht)

    # hmm TODO
    initial_size = 2^3

    exponents = Vector{Vector{UInt16}}(undef, initial_size)
    hashdata  = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(Int, initial_size)

    # preserve ring info
    explen = basis_ht.explen
    nvars  = basis_ht.nvars
    ord    = basis_ht.ord

    # preserve division info
    divmap   = basis_ht.divmap
    divvars  = basis_ht.divvars
    ndivvars = basis_ht.ndivvars
    ndivbits = basis_ht.ndivbits

    # TODO: preserve hasher?
    hasher = basis_ht.hasher

    load = 1
    size = initial_size
    offset = 2

    exponents[1] = zeros(UInt16, explen)

    CustomMonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

#------------------------------------------------------------------------------

# resizes (if needed) ht so that it can store `size` elements,
# and clears all previoud data
function reinitialize_hash_table!(ht::CustomMonomialHashtable, size::Int)
    if size > ht.size
        while size > ht.size
            ht.size *= 2
        end
        # TODO
        resize!(ht.hashdata, ht.size)
        resize!(ht.exponents, ht.size)
    end
    ht.hashtable = zeros(Int, ht.size)
    ht.load = 1
end

# doubles the size of storage in `ht`,
# and rehashes all elements
function enlarge_hash_table!(ht::CustomMonomialHashtable)
    ht.size *= 2
    resize!(ht.hashdata,  ht.size)
    resize!(ht.exponents, ht.size)

    ht.hashtable = zeros(Int, ht.size)

    mod = ht.size - 1
    for i in ht.offset:ht.load  # TODO: hardcoding 2 is bad
        # hash for this elem is already computed
        he = ht.hashdata[i].hash
        hidx = he
        for j in 0:ht.size-1
            hidx = (hidx + j) % mod + 1
            ht.hashtable[hidx] != 0 && continue
            ht.hashtable[hidx] = i
            break
        end
    end
end

#------------------------------------------------------------------------------

function insert_in_hash_table!(ht::CustomMonomialHashtable, e::Vector{UInt16})
    # generate hash
    he = UInt32(0)
    for i in 1:ht.explen
        he += ht.hasher[i] * e[i]
    end

    @debug "new hash:" he

    # find new elem position in the table
    hidx = he
    mod = ht.size

    for i in 0:ht.size - 1
        hidx = (hidx + i) % mod + 1
        vidx  = ht.hashtable[hidx]

        @info "" i hidx vidx

        # if free
        vidx == 0 && break

        # if not free and not same hash
        if ht.hashdata[vidx].hash != he
            continue
        end

        present = ht.exponents[vidx]
        for j in 1:ht.explen
            # if hash collision
            if present[j] != e[j]
                continue
            end
        end

        # already present in hashtable
        return vidx
    end

    # add its position to hashtable, and insert exponent to that position
    vidx = ht.load + 1
    ht.hashtable[hidx] = vidx
    ht.exponents[vidx] = e

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
function fill_divmask!(ht::CustomMonomialHashtable)
    ndivvars = ht.ndivvars
    divvars  = ht.divvars

    min_exp  = Vector{UInt16}(undef, ndivvars)
    max_exp  = Vector{UInt16}(undef, ndivvars)

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

    @info "filling divmask" min_exp max_exp ht.ndivbits

    ctr   = 1
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

    @info "filling divmask" ht.divmap

    # TODO: offset
    for vidx in ht.offset:ht.load
        unmasked = ht.hashdata[vidx]
        e = ht.exponents[vidx]
        divmask  = generate_monomial_divmask(e, ht)
        ht.hashdata[vidx] = Hashvalue(unmasked.hash, divmask, 0, e[end])
    end
end

#=
    TODO

=#
function generate_monomial_divmask(
            e::Vector{UInt16},
            ht::CustomMonomialHashtable)

    divvars = ht.divvars
    divmap  = ht.divmap

    # TODO: we miss one bit of space here
    # upd: not anymore
    ctr = UInt32(1)
    res = UInt32(0)
    for i in 1:ht.ndivvars
        for j in 1:ht.ndivbits
            if e[divvars[i]] >= divmap[ctr]
                res |= UInt32(1) << (ctr - 1)
            end
            ctr += 1
        end
    end

    res
end

#------------------------------------------------------------------------------

# h1 divisible by h2
function is_monom_divisible(h1::Int, h2::Int, ht::CustomMonomialHashtable)

    if (ht.hashdata[h2].divmask & ~ht.hashdata[h1].divmask) != 0
        @debug "used divmask"
        return false
    end

    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]
    # TODO: one less iteration is possible
    for i in 1:ht.explen
        if e1[i] < e2[i]
            return false
        end
    end

    return true
end

function is_gcd_const(h1::Int, h2::Int, ht::CustomMonomialHashtable)
    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]

    # TODO: explen - 1
    for i in 1:ht.explen - 1
        if e1[i] != 0 && e2[i] != 0
            return false
        end
    end

    return true
end

#------------------------------------------------------------------------------

function check_monomial_division_in_update(a, first, last, lcm, ht)
    # pairs are sorted, we only need to check entries above starting point

    divmask = ht.hashdata[lcm].divmask
    lcmexp  = ht.exponents[lcm]

    j = first + 1
    @label Restart
    while j < last
        a[j] == 0 && continue
        (~ht.hashdata[a[j]].divmask & divmask != 0) && continue
        ea = ht.exponents[a[j]]
        for i in 1:ht.explen
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
# with hash `htmp` to hashtable `symbol_ht`
# Who is `row` then?
function insert_multiplied_poly_in_hash_table!(
        row,
        htmp,
        etmp,
        poly,
        ht,
        symbol_ht)

    # oof

    # length of poly to add
    len    = length(poly)
    explen = ht.explen

    mod = ht.size - 1

    bexps = ht.exponents
    bdata = ht.hashdata

    sexps = symbol_ht.exponents
    sdata = symbol_ht.hashdata

    @info "adding poly $poly mult by $etmp to symbolic hashtable"

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
        h = htmp + bdata[poly[l]].hash
        e = bexps[poly[l]]
        println("monom of index $(poly[l]) : $e")

        lastidx = symbol_ht.load + 1
        sexps[lastidx] = Vector{UInt16}(undef, explen)
        enew = sexps[lastidx]
        for j in 1:explen
            # multiplied monom exponent
            enew[j] = etmp[j] + e[j]
        end

        @info "product is $enew"
        # now insert into hashtable
        k = h
        i = 0
        @label Restart
        while i <= symbol_ht.size  # TODO: < or <= ?
            k = (k + i) % mod + 1
            vidx = symbol_ht.hashtable[k]
            # if index is free
            vidx == 0 && break
            # if different exponent is stored here
            sdata[vidx].hash != h && continue

            estored = sexps[vidx]
            for j in 1:explen
                # hash collision, restarting search
                if estored[j] != enew[j]
                    i += 1
                    @goto Restart
                end
            end

            row[l] = vidx
            l += 1

            @goto Letsgo
        end

        # add multiplied exponent to hash table
        symbol_ht.hashtable[k] = lastidx

        divmask = generate_monomial_divmask(enew, symbol_ht)
        sdata[lastidx] = Hashvalue(h, divmask, 0, enew[end])

        row[l] = lastidx

        symbol_ht.load += 1
    end

    row
end

function multiplied_poly_to_matrix_row!(symbolic_ht, basis_ht, htmp, etmp, poly)
    row = copy(poly)
    while symbolic_ht.load + length(poly) >= symbolic_ht.size
        enlarge_hash_table!(symbolic_ht)
    end
    insert_multiplied_poly_in_hash_table!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
end

#------------------------------------------------------------------------------

function insert_in_basis_hash_table_pivots(
        row,
        ht,
        symbol_ht,
        col2hash)

    while ht.size - ht.load <= length(row)
        enlarge_hash_table!(ht)
    end

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.exponents

    mod    = ht.size - 1
    explen = ht.explen
    bdata  = ht.hashdata
    bexps  = ht.exponents
    bhash  = ht.hashtable

    l = 1
    @label Letsgo
    while l <= length(row)
        hidx = col2hash[row[l]]
        h    = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = copy(sexps[hidx])
        e = bexps[lastidx]

        k = h
        i = 0
        @label Restart
        while i < ht.size
            k = (k + i) % mod + 1
            hm = bhash[k]

            hm == 0 && break
            bdata[hm].hash != h && continue

            ehm = exps[hm]
            for j in 1:explen
                if e[j] != ehm[j]
                    i += 1
                    @goto Restart
                end
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        hmap[k] = pos = lastidx
        row[l] = pos

        bdata[k] = Hashvalue(h, sdata[hidx].divmask,
                             sdata[hidx].idx, sdata[hidx].deg)

        ht.load += 1
    end

end

function insert_plcms_in_basis_hash_table!(
        pairset,
        off,
        ht,
        update_ht,
        basis,
        plcm,
        ifirst, ilast)

    gens  = basis.gens
    mod   = ht.size
    ps    = pairset.pairs

    @info "insert plcms in basis hash table" plcm ifirst ilast

    m = ifirst
    l = 1
    @label Letsgo
    while l <= ilast
        # what? why??
        if is_gcd_const(gens[ps[off + l].poly1][1], gens[ps[off + 1].poly2][1], ht)
            continue
        end

        ps[m] = ps[off + l]

        @debug "hmm" ht.load ht.size plcm[l] update_ht.size update_ht.load

        h = update_ht.hashdata[plcm[l]].hash
        ht.exponents[ht.load+1] = copy(update_ht.exponents[plcm[l]])
        n = ht.exponents[ht.load+1]

        k = h
        i = 0
        @label Restart
        while i < ht.size
            k = (k + i) % mod
            hm = ht.hashtable[k]

            hm != 0 && break
            ht.hashdata[hm].hash != h && continue

            ehm = ht.exponents[hm]
            for j in 1:ht.explen
                if ehm[j] != n[j]
                    i += 1
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
        ht.hashdata[ht.load + 1] = Hashvalue(h, uhd[ll].divmask, 0, uhd[ll].deg)

        ht.load += 1
        ps[m] = SPair(ps[m].poly1, ps[m].poly2, pos, ps[m].deg)
        m += 1
    end

end

#------------------------------------------------------------------------------

# computes lcm of he1 and he2 and inserts it in ht2
function get_lcm(he1, he2, ht1, ht2)
    e1   = ht1.exponents[he1]
    e2   = ht1.exponents[he2]
    etmp = ht1.exponents[1]

    # TODO: degrevlex only
    etmp[end] = 0
    for i in 1:ht1.explen-1
        etmp[i]    = max(e1[i], e2[i])
        etmp[end] += etmp[i]
    end

    insert_in_hash_table!(ht2, etmp)
end

#------------------------------------------------------------------------------
