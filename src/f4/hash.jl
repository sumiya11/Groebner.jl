

# TODO: make all insertions to be carried via one method (?)
#       Is that possible? Yes

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

#=
    Hashtable implementing linear probing
    and designed to store and operate with monomials
=#
mutable struct MonomialHashtable
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
function initialize_basis_hash_table(
        ring::PolyRing,
        rng::Random.AbstractRNG;
        initial_size::Int=2^16)

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
            # TODO: pass random generator
            hasher[i] = rand(rng, UInt32)
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
    charbit    = 8 # TODO ??
    int32bits  = charbit * sizeof(Int32)
    int32bits != 32 && error("Strange story with ints")
    ndivbits   = div(int32bits, nvars)
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

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

# initialize hashtable either for `symbolic_preprocessing` or for `update` functions
# These are of the same purpose and structure as basis hashtable,
# but are more local oriented
function initialize_secondary_hash_table(basis_ht::MonomialHashtable)

    # hmm TODO
    initial_size = 2^6

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
    # yes
    hasher = basis_ht.hasher

    load = 1
    size = initial_size
    offset = 2

    exponents[1] = zeros(UInt16, explen)

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

#------------------------------------------------------------------------------

# resizes (if needed) ht so that it can store `size` elements,
# and clears all previoud data
function reinitialize_hash_table!(ht::MonomialHashtable, size::Int)
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
function enlarge_hash_table!(ht::MonomialHashtable)
    ht.size *= 2
    resize!(ht.hashdata,  ht.size)
    resize!(ht.exponents, ht.size)

    ht.hashtable = zeros(Int, ht.size)

    mod = UInt32(ht.size - 1)
    for i in ht.offset:ht.load  # TODO: hardcoding 2 is bad
        # hash for this elem is already computed
        he = ht.hashdata[i].hash
        hidx = he
        for j in UInt32(1):UInt32(ht.size)
            hidx = (he + j) & mod + UInt32(1)
            ht.hashtable[hidx] != 0 && continue
            ht.hashtable[hidx] = i
            break
        end
    end
end

#------------------------------------------------------------------------------

function insert_in_hash_table!(ht::MonomialHashtable, e::Vector{UInt16})
    # generate hash
    he = UInt32(0)
    for i in 1:ht.explen
        he += ht.hasher[i] * e[i]
    end
    # find new elem position in the table
    hidx = Int(he)  # Int for type stability
    # power of twoooo
    mod = UInt32(ht.size - 1)
    i = UInt32(1)

    @label Restart
    while i < ht.size
        # TODO: & instead of %
        hidx = (he + i) & mod + UInt32(1)
        @inbounds vidx  = ht.hashtable[hidx]

        # if free
        vidx == 0 && break

        # if not free and not same hash
        if ht.hashdata[vidx].hash != he
            i += UInt32(1)
            continue
        end

        present = ht.exponents[vidx]
        for j in 1:ht.explen
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
    ht.exponents[vidx] = copy(e)

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
            ht::MonomialHashtable)

    divvars = ht.divvars
    divmap  = ht.divmap

    ctr = UInt32(1)
    res = UInt32(0)
    for i in 1:ht.ndivvars
        for j in 1:ht.ndivbits
            @inbounds if e[divvars[i]] >= divmap[ctr]
                res |= UInt32(1) << (ctr - 1)
            end
            ctr += UInt32(1)  # for typ stability
        end
    end

    res
end

#------------------------------------------------------------------------------

# h1 divisible by h2
function is_monom_divisible(h1::Int, h2::Int, ht::MonomialHashtable)

    if (ht.hashdata[h2].divmask & ~ht.hashdata[h1].divmask) != 0
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

function is_gcd_const(h1::Int, h2::Int, ht::MonomialHashtable)
    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]

    for i in 1:ht.explen - 1
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
    lcmexp  = ht.exponents[lcm]

    j = first
    @label Restart
    while j <= last
        # bad lcm
        if a[j] == 0
            j += 1
            continue
        end
        # fast division check
        if (~ht.hashdata[a[j]].divmask & divmask) != 0
            j += 1
            continue
        end
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
# with hash `htmp` to hashtable `symbol_ht`,
# and substitute hashes in row
function insert_multiplied_poly_in_hash_table!(
        row::Vector{Int},
        htmp::UInt32,
        etmp::Vector{UInt16},
        poly::Vector{Int},
        ht::MonomialHashtable,
        symbol_ht::MonomialHashtable)

    # oof

    # length of poly to add
    len    = length(poly)
    explen = ht.explen

    mod = UInt32(symbol_ht.size - 1)

    bexps = ht.exponents
    bdata = ht.hashdata

    sexps = symbol_ht.exponents
    sdata = symbol_ht.hashdata

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
        # println("monom of index $(poly[l]) : $e")

        lastidx = symbol_ht.load + 1
        if !isassigned(sexps, lastidx)
            sexps[lastidx] = Vector{UInt16}(undef, explen)
        end
        enew = sexps[lastidx]
        @inbounds for j in 1:explen
            # multiplied monom exponent
            enew[j] = etmp[j] + e[j]
        end

        # now insert into hashtable
        k = h
        i = UInt32(1)
        @label Restart
        while i <= symbol_ht.size  # TODO: < or <= ?
            # @assert ispow2(mod + 1)
            k = (h + i) & mod + UInt32(1)
            @inbounds vidx = symbol_ht.hashtable[k]
            # if index is free
            vidx == 0 && break
            # if different exponent is stored here
            if sdata[vidx].hash != h
                i += UInt32(1)
                continue
            end

            estored = sexps[vidx]
            @inbounds for j in 1:explen
                # hash collision, restarting search
                if estored[j] != enew[j]
                    i += UInt32(1)
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
        l += 1
        symbol_ht.load += 1
    end

    row
end

function multiplied_poly_to_matrix_row!(
        symbolic_ht::MonomialHashtable, basis_ht::MonomialHashtable,
        htmp::UInt32, etmp::Vector{UInt16}, poly::Vector{Int})

    row = copy(poly)
    while symbolic_ht.load + length(poly) >= symbolic_ht.size
        enlarge_hash_table!(symbolic_ht)
    end
    insert_multiplied_poly_in_hash_table!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
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

    mod    = UInt32(ht.size - 1)
    explen = ht.explen
    bdata  = ht.hashdata
    bexps  = ht.exponents
    bhash  = ht.hashtable

    l = 1
    @label Letsgo
    while l <= length(row)
        hidx = col2hash[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        # TODO: speed this up
        bexps[lastidx] = copy(sexps[hidx])
        e = bexps[lastidx]

        k = h
        i = UInt32(1)
        @label Restart
        while i <= ht.size
            k = (h + i) & mod + UInt32(1)
            hm = bhash[k]

            hm == 0 && break
            if bdata[hm].hash != h
                i += UInt32(1)
                continue
            end

            ehm = bexps[hm]
            for j in 1:explen
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

    gens  = basis.gens
    mod   = UInt32(ht.size - 1)
    ps    = pairset.pairs

    m = ifirst
    l = 1
    @label Letsgo
    while l < ilast
        if plcm[l] == 0
            l += 1
            continue
        end

        if is_gcd_const(gens[ps[off + l].poly1][1], gens[ps[off + 1].poly2][1], ht)
            l += 1
            continue
        end

        ps[m] = ps[off + l]

        # TODO: IT IS NOT CORRECT
        # upd: it is, but it can be done better
        h = update_ht.hashdata[plcm[l]].hash
        ht.exponents[ht.load+1] = copy(update_ht.exponents[plcm[l]])
        n = ht.exponents[ht.load+1]

        k = h
        i = UInt32(1)
        @label Restart
        while i <= ht.size
            k = (h + i) & mod + UInt32(1)
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
        l += 1
    end

    pairset.load = m - 1
end

#------------------------------------------------------------------------------

# computes lcm of he1 and he2 as exponent vectors from ht1
# and inserts it in ht2
function get_lcm(he1::Int, he2::Int,
                    ht1::MonomialHashtable, ht2::MonomialHashtable)

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
