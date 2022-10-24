
#------------------------------------------------------------------------------
####### Hashtable ######

#=
    stores hash of a monomial,
    corresponding divmask to speed up divisibility checks,
    index for position matrix (defaults to zero),
    and the todal degree
=#
mutable struct Hashvalue
    hash::ExponentHash
    divmask::DivisionMask
    idx::Int
    deg::Degree
end

function copy_hashvalue(x::Hashvalue)
    Hashvalue(x.hash, x.divmask, x.idx, x.deg)
end

#=
    Hashtable implementing linear probing
    and designed to store and operate with monomials
=#
mutable struct MonomialHashtable
    exponents::Vector{ExponentVector}

    # maps exponent hash to its position in exponents array
    hashtable::Vector{ExponentIdx}

    # stores hashes, division masks,
    # and other valuable info
    # for each hashtable enrty
    hashdata::Vector{Hashvalue}

    # values to hash exponents with, i.e
    # hash(e) = sum(hasher .* e)
    hasher::Vector{ExponentHash}

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

# Returns the next look-up position in the table 
# (that is, implementing open addressing with quadratic probing) 
function hashnextindex(h::ExponentHash, j::ExponentHash, mod::ExponentHash)
    (h + j) & mod + ExponentHash(1)
end

#------------------------------------------------------------------------------

# initialize and set fields for basis hashtable
function initialize_basis_hash_table(
    ring::PolyRing,
    rng::Random.AbstractRNG;
    initial_size::Int=2^16)

    # not necessary to create `initial_size` exponents
    exponents = Vector{ExponentVector}(undef, initial_size)
    hashdata = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(ExponentIdx, initial_size)

    nvars = ring.nvars
    explen = ring.explen
    ord = ring.ord

    # initialize hashing vector
    hasher = zeros(ExponentHash, explen)
    for i in 1:explen
        # we don't want hash vector components to be zero
        while iszero(hasher[i])
            hasher[i] = rand(rng, ExponentHash)
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
    charbit = 8
    int32bits = charbit * sizeof(Int32)
    int32bits != 32 && error("Strange story with ints")
    ndivbits = div(int32bits, nvars)
    # division mask stores at least 1 bit
    # per each of first charbit*sizeof(Int32) variables
    ndivbits == 0 && (ndivbits += 1)
    ndivvars = nvars < int32bits ? nvars : int32bits
    divvars = Vector{Int}(undef, ndivvars)
    divmap = Vector{UInt32}(undef, ndivvars * ndivbits)
    # count only first ndivvars variables for divisibility checks
    for i in 1:ndivvars
        divvars[i] = i
    end

    println("On HT initialize");
    println(ndivvars, " ", ndivbits);
    println(divvars, " ", divmap);

    # first stored exponent used as buffer lately
    exponents[1] = zeros(Degree, explen)

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, explen, ord,
        divmap, divvars, ndivvars, ndivbits,
        size, load, offset)
end

function copy_hashtable(ht::MonomialHashtable)
    exps = Vector{ExponentVector}(undef, ht.size)
    table = Vector{ExponentIdx}(undef, ht.size)
    data = Vector{Hashvalue}(undef, ht.size)
    exps[1] = zeros(Degree, ht.explen)

    @inbounds for i in 2:ht.load
        exps[i] = copy(ht.exponents[i])
        table[i] = ht.hashtable[i]
        data[i] = copy_hashvalue(ht.hashdata[i])
    end

    MonomialHashtable(
        ht.exponents, ht.hashtable, ht.hashdata, ht.hasher,
        ht.nvars, ht.explen, ht.ord,
        ht.divmap, ht.divvars, ht.ndivvars, ht.ndivbits,
        ht.size, ht.load, ht.offset)
end

#------------------------------------------------------------------------------

# initialize hashtable either for `symbolic_preprocessing` or for `update` functions
# These are of the same purpose and structure as basis hashtable,
# but are more local oriented
function initialize_secondary_hash_table(basis_ht::MonomialHashtable)

    # 2^6 seems to be the best out of 2^5, 2^6, 2^7
    initial_size = 2^6

    exponents = Vector{ExponentVector}(undef, initial_size)
    hashdata = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(ExponentIdx, initial_size)

    # preserve ring info
    explen = basis_ht.explen
    nvars = basis_ht.nvars
    ord = basis_ht.ord

    # preserve division info
    divmap = basis_ht.divmap
    divvars = basis_ht.divvars
    ndivvars = basis_ht.ndivvars
    ndivbits = basis_ht.ndivbits

    # preserve hasher
    hasher = basis_ht.hasher

    load = 1
    size = initial_size
    offset = 2

    exponents[1] = zeros(Degree, explen)

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
        resize!(ht.hashdata, ht.size)
        resize!(ht.exponents, ht.size)
    end
    ht.hashtable = zeros(Int, ht.size)
    ht.load = 1
end

# doubles the size of storage in `ht`,
# and rehashes all elements
function enlarge_hash_table!(ht::MonomialHashtable)
    # probably not correct

    ht.size *= 2
    resize!(ht.hashdata, ht.size)
    resize!(ht.exponents, ht.size)

    ht.hashtable = zeros(Int, ht.size)

    mod = ExponentHash(ht.size - 1)
    for i in ht.offset:ht.load
        # hash for this elem is already computed
        he = ht.hashdata[i].hash
        hidx = he
        @inbounds for j in UInt32(1):UInt32(ht.size)
            hidx = hashnextindex(he, j, mod)
            ht.hashtable[hidx] != 0 && continue
            ht.hashtable[hidx] = i
            break
        end
    end
end


#------------------------------------------------------------------------------

function insert_in_hash_table!(ht::MonomialHashtable, e::ExponentVector)
    # generate hash
    he = ExponentHash(0)

    # here, e[i] is of type UInt16, while hasher[i] is UInt32 =(
    @inbounds for i in 1:ht.explen
        he += ht.hasher[i] * e[i]
    end

    # find new elem position in the table
    hidx = ExponentHash(he)
    # power of twoooo
    @assert ispow2(ht.size)
    mod = ExponentHash(ht.size - 1)
    i = ExponentHash(1)

    @label Restart
    while i < ht.size
        hidx = hashnextindex(he, i, mod)
        @inbounds vidx = ht.hashtable[hidx]

        # if free
        iszero(vidx) && break

        # if not free and not same hash
        if ht.hashdata[vidx].hash != he
            i += ExponentHash(1)
            continue
        end

        present = ht.exponents[vidx]
        @inbounds for j in 1:ht.explen
            # if hash collision
            if present[j] != e[j]
                i += ExponentHash(1)
                @goto Restart
            end
        end

        # already present in hashtable
        return vidx
    end

    # add its position to hashtable, and insert exponent to that position
    vidx = ExponentIdx(ht.load + 1)
    ht.hashtable[hidx] = vidx

    # probably can be changed 

    # TODO: check efficiency
    ht.exponents[vidx] = similar(e)
    ve = ht.exponents[vidx]
    @inbounds for i in 1:length(e)
        ve[i] = e[i]
    end

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

    min_exp = ExponentVector(undef, ndivvars)
    max_exp = ExponentVector(undef, ndivvars)

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

    println("In fill_divmask")
    println(ht.ndivbits, " ", Int.(min_exp), " ", Int.(max_exp))

    ctr = 1
    steps = UInt32(0)
    for i in 1:ndivvars
        steps = div(max_exp[i] - min_exp[i], UInt32(ht.ndivbits))
        (iszero(steps)) && (steps += UInt32(1))
        for j in 1:ht.ndivbits
            ht.divmap[ctr] = steps
            steps += UInt32(1)
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
    e::ExponentVector,
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
function is_monom_divisible(h1::ExponentIdx, h2::ExponentIdx, ht::MonomialHashtable)
    @inbounds if !iszero(ht.hashdata[h2].divmask & ~ht.hashdata[h1].divmask)
        return false
    end

    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]
    @inbounds for i in 1:ht.explen
        if e1[i] < e2[i]
            return false
        end
    end

    return true
end

function is_gcd_const(h1::ExponentIdx, h2::ExponentIdx, ht::MonomialHashtable)
    e1 = ht.exponents[h1]
    e2 = ht.exponents[h2]

    for i in 1:ht.explen-1
        if !iszero(e1[i]) && !iszero(e2[i])
            return false
        end
    end

    return true
end

#------------------------------------------------------------------------------

# compare pairwise divisibility of lcms from a[first:last] with lcm
function check_monomial_division_in_update(
    a::Vector{ExponentIdx}, first::Int, last::Int,
    lcm::ExponentIdx, ht::MonomialHashtable)

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
    row::Vector{ExponentIdx},
    htmp::ExponentHash,
    etmp::ExponentVector,
    poly::Vector{ExponentIdx},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable)

    # oof

    # length of poly to add
    len = length(poly)
    explen = ht.explen

    mod = ExponentHash(symbol_ht.size - 1)

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

        i = ExponentHash(1)

        @label Restart
        while i <= symbol_ht.size
            k = hashnextindex(h, i, mod)

            @inbounds vidx = symbol_ht.hashtable[k]
            # if index is free
            iszero(vidx) && break
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
                    i += ExponentHash(1)

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

# If symbolic hash table of the given size and load factor should be
# enlarged after adding the polynomial of length `added` 
function symbolic_ht_needscale(load::Int, added::Int, size::Int)
    1.4*(load + added) >= size
end

function multiplied_poly_to_matrix_row!(
    symbolic_ht::MonomialHashtable, basis_ht::MonomialHashtable,
    htmp::ExponentHash, etmp::ExponentVector, poly::MonomsVector)

    row = similar(poly)
    while symbolic_ht_needscale(symbolic_ht.load, length(poly), symbolic_ht.size)
        enlarge_hash_table!(symbolic_ht)
    end

    insert_multiplied_poly_in_hash_table!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
end

#------------------------------------------------------------------------------

function insert_in_basis_hash_table_pivots(
    row::Vector{ColumnIdx},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    col2hash::Vector{ExponentIdx})

    while ht.size - ht.load <= length(row)
        enlarge_hash_table!(ht)
    end

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.exponents

    mod = ExponentHash(ht.size - 1)
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
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = ExponentHash(1)
        @label Restart
        while i <= ht.size
            k = hashnextindex(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break
            @inbounds if bdata[hm].hash != h
                i += ExponentHash(1)
                continue
            end

            ehm = bexps[hm]
            @inbounds for j in 1:explen
                if e[j] != ehm[j]
                    i += ExponentHash(1)
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

#------------------------------------------------------------------------------

# computes lcm of he1 and he2 as exponent vectors from ht1
# and inserts it in ht2
function get_lcm(he1::ExponentIdx, he2::ExponentIdx,
    ht1::MonomialHashtable, ht2::MonomialHashtable)

    @inbounds e1 = ht1.exponents[he1]
    @inbounds e2 = ht1.exponents[he2]
    @inbounds etmp = ht1.exponents[1]

    @inbounds etmp[end] = 0
    @inbounds for i in 1:ht1.explen-1
        etmp[i] = max(e1[i], e2[i])
        etmp[end] += etmp[i]
    end

    insert_in_hash_table!(ht2, etmp)
end