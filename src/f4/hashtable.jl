# Hashtable

# This hashtable stores monomials.
# This hashtable implementation assumes 
# that the hash function that acts on monomials is linear. 
# (each monomial type should implement linear hash function)

# Some monomial implementations are mutable, and some are not.
# In order to maintain generic enough code that will work for both
# we usually write something like:
#   enew = monom_product!(enew, e1, e2)
# Then, if a mutable implementation is used, 
# enew would be overwritten inside of monom_product!.
# Otherwise, monom_product! would return a new immutable object
# that is then assigned to enew.
# That allows us to write more or less independently
# of the monomial implementation.

# The hashtable size is always the power of two.
# The hashtable size is doubled each time the load factor exceeds
#   ht_resize_threshold()
# (this ht_resize_threshold() should be around 0.5
#  to balance the rehashing cost with the insertion cost)
# (this ht_resize_threshold() should be a bit smaller than 0.5,
#  since the number of hit insertions greatly exceeds the number
#  of miss insertions in the hashtable)

# Hashvalue.
# stores index for position in the matrix (defaults to zero),
# hash of a monomial,
# corresponding divmask to speed up divisibility checks,
# and the todal degree
mutable struct Hashvalue
    idx::Int
    hash::MonomHash
    divmask::DivisionMask
    deg::MonomHash
end

function copy_hashvalue(x::Hashvalue)
    Hashvalue(x.idx, x.hash, x.divmask, x.deg)
end

# Hashtable designed to store monomials
mutable struct MonomialHashtable{M<:Monom, Ord<:AbstractMonomialOrdering}
    exponents::Vector{M}

    # maps exponent hash to its position in exponents array
    hashtable::Vector{MonomIdx}

    # stores hashes, division masks,
    # and other valuable info
    # for each hashtable enrty
    hashdata::Vector{Hashvalue}

    # values to hash exponents with, i.e
    # hash(e) := hash(e, hasher)
    hasher::Vector{MonomHash}

    #= Ring information =#
    # number of variables
    nvars::Int
    # ring monomial ordering
    ord::Ord

    #= Monom divisibility =#
    # divisor map to check divisibility faster
    divmap::Vector{UInt32}
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
function nexthashindex(h::MonomHash, j::MonomHash, mod::MonomHash)
    (h + j) & mod + MonomHash(1)
end

#------------------------------------------------------------------------------

# initialize and set fields for basis hashtable
function initialize_basis_hash_table(
    ring::PolyRing{Char, Ord},
    rng::Random.AbstractRNG,
    MonomT;
    initial_size::Int=2^16) where {Char, Ord<:AbstractMonomialOrdering}

    # not necessary to create `initial_size` exponents
    exponents = Vector{MonomT}(undef, initial_size)
    hashdata = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(MonomIdx, initial_size)

    nvars = ring.nvars
    ord = ring.ord

    # initialize hashing vector
    hasher = make_hasher(MonomT, nvars)

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
    int32bits != 32 && error("Strange story with Ints")
    ndivbits = div(int32bits, nvars)
    # division mask stores at least 1 bit
    # per each of first ndivvars variables
    ndivbits == 0 && (ndivbits += 1)
    # count only first ndivvars variables for divisibility checks
    ndivvars = nvars < int32bits ? nvars : int32bits
    divmap = Vector{DivisionMask}(undef, ndivvars * ndivbits)

    # first stored exponent used as buffer lately
    exponents[1] = make_zero_ev(MonomT, nvars)

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, ord,
        divmap, ndivvars, ndivbits,
        size, load, offset)
end

function copy_hashtable(ht::MonomialHashtable{M, O}) where {M, O}
    exps = Vector{M}(undef, ht.size)
    table = Vector{MonomIdx}(undef, ht.size)
    data = Vector{Hashvalue}(undef, ht.size)
    exps[1] = make_zero_ev(M, ht.nvars)

    @inbounds for i in 2:ht.load
        exps[i] = copy(ht.exponents[i])
        table[i] = ht.hashtable[i]
        data[i] = copy_hashvalue(ht.hashdata[i])
    end

    MonomialHashtable(
        ht.exponents, ht.hashtable, ht.hashdata, ht.hasher,
        ht.nvars, ht.ord,
        ht.divmap, ht.ndivvars, ht.ndivbits,
        ht.size, ht.load, ht.offset)
end

# initialize hashtable either for `symbolic_preprocessing` or for `update` functions
# These are of the same purpose and structure as basis hashtable,
# but are more local oriented
function initialize_secondary_hash_table(basis_ht::MonomialHashtable{M}) where {M}
    # 2^6 seems to be the best out of 2^5, 2^6, 2^7
    initial_size = 2^6

    exponents = Vector{M}(undef, initial_size)
    hashdata = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(MonomIdx, initial_size)

    # preserve ring info
    nvars = basis_ht.nvars
    ord = basis_ht.ord

    # preserve division info
    divmap = basis_ht.divmap
    ndivbits = basis_ht.ndivbits
    ndivvars = basis_ht.ndivvars

    # preserve hasher
    hasher = basis_ht.hasher

    load = 1
    size = initial_size
    offset = 2

    exponents[1] = make_zero_ev(M, nvars)

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, ord,
        divmap, ndivvars, ndivbits,
        size, load, offset)
end

function select_tablesize(ring::PolyRing, polys::AbstractVector)
    nvars = ring.nvars
    sz = length(polys)

    tablesize = 2^10
    if nvars > 4
        tablesize = 2^14
    end
    if nvars > 7
        tablesize = 2^16
    end

    if sz < 3
        tablesize = div(tablesize, 2)
    end
    if sz < 2
        tablesize = div(tablesize, 2)
    end

    tablesize
end

#------------------------------------------------------------------------------

ht_resize_threshold() = 0.4
ht_needs_resize(size, load, added) = (load + added)/size > ht_resize_threshold()

function check_enlarge_hashtable!(ht::MonomialHashtable, added::Integer)
    newsize = ht.size
    while ht_needs_resize(newsize, ht.load, added)
        newsize *= 2
    end
    if newsize != ht.size
        ht.size = newsize
        @assert ispow2(ht.size)

        resize!(ht.hashdata, ht.size)
        resize!(ht.exponents, ht.size)
        ht.hashtable = zeros(Int, ht.size)
        
        mod = MonomHash(ht.size - 1)

        for i in ht.offset:ht.load
            # hash for this elem is already computed
            he = ht.hashdata[i].hash
            hidx = he
            @inbounds for j in MonomHash(1):MonomHash(ht.size)
                hidx = nexthashindex(he, j, mod)
                !iszero(ht.hashtable[hidx]) && continue
                ht.hashtable[hidx] = i
                break
            end
        end
    end
    nothing
end

#------------------------------------------------------------------------------

# if hash collision happened
function ishashcollision(ht::MonomialHashtable, vidx, e, he)
    # if not free and not same hash
    @inbounds if ht.hashdata[vidx].hash != he
        return true
    end
    # if not free and not same monomial
    @inbounds if !is_monom_elementwise_eq(ht.exponents[vidx], e)
        return true
    end
    false
end

function insert_in_hash_table!(ht::MonomialHashtable{M}, e::M) where {M}
    # generate hash
    he::MonomHash = hash(e, ht.hasher)

    # find new elem position in the table
    hidx = MonomHash(he)
    # power of twoooo
    @assert ispow2(ht.size)
    mod = MonomHash(ht.size - 1)
    i = MonomHash(1)
    hsize = MonomHash(ht.size)

    @inbounds while i < hsize
        hidx = nexthashindex(he, i, mod)

        vidx = ht.hashtable[hidx]

        # if free
        iszero(vidx) && break

        # if not free and not same hash
        if ishashcollision(ht, vidx, e, he)
            i += MonomHash(1)
            continue
        end

        # already present in hashtable
        return vidx
    end

    # add its position to hashtable, and insert exponent to that position
    vidx = MonomIdx(ht.load + 1)
    @inbounds ht.hashtable[hidx] = vidx
    @inbounds ht.exponents[vidx] = copy(e)
    divmask = monom_divmask(e, DivisionMask, ht.ndivvars, ht.divmap, ht.ndivbits)
    @inbounds ht.hashdata[vidx] = Hashvalue(0, he, divmask, totaldeg(e))

    ht.load += 1

    return vidx
end

#------------------------------------------------------------------------------

function is_divmask_divisible(d1::DivisionMask, d2::DivisionMask)
    iszero(~d1 & d2)
end

#=
    Having `ht` filled with monomials from input polys,
    computes ht.divmap and divmask for each of already stored monomials
=#
function fill_divmask!(ht::MonomialHashtable)
    ndivvars = ht.ndivvars

    min_exp = Vector{UInt64}(undef, ndivvars)
    max_exp = Vector{UInt64}(undef, ndivvars)

    e = Vector{UInt64}(undef, ht.nvars)
    make_dense!(e, ht.exponents[ht.offset])

    @inbounds for i in 1:ndivvars
        min_exp[i] = e[i]
        max_exp[i] = e[i]
    end

    @inbounds for i in ht.offset:ht.load # TODO: offset
        make_dense!(e, ht.exponents[i])
        for j in 1:ndivvars
            if e[j] > max_exp[j]
                max_exp[j] = e[j]
                continue
            end
            if e[j] < min_exp[j]
                min_exp[j] = e[j]
            end
        end
    end

    ctr = 1
    steps = UInt32(0)
    @inbounds for i in 1:ndivvars
        steps = div(max_exp[i] - min_exp[i], UInt32(ht.ndivbits))
        (iszero(steps)) && (steps += UInt32(1))
        for j in 1:ht.ndivbits
            ht.divmap[ctr] = steps
            steps += UInt32(1)
            ctr += 1
        end
    end
    @inbounds for vidx in ht.offset:ht.load
        unmasked = ht.hashdata[vidx]
        e = ht.exponents[vidx]
        divmask = monom_divmask(e, DivisionMask, ht.ndivvars, ht.divmap, ht.ndivbits)
        ht.hashdata[vidx] = Hashvalue(0, unmasked.hash, divmask, totaldeg(e))
    end

    nothing
end

#------------------------------------------------------------------------------

# h1 divisible by h2
function is_monom_divisible(h1::MonomIdx, h2::MonomIdx, ht::MonomialHashtable)
    @inbounds if !is_divmask_divisible(ht.hashdata[h1].divmask, ht.hashdata[h2].divmask)
        return false
    end
    @inbounds e1 = ht.exponents[h1]
    @inbounds e2 = ht.exponents[h2]
    is_monom_divisible(e1, e2)
end

# checks that gcd(g1, h2) is one
function is_gcd_const(h1::MonomIdx, h2::MonomIdx, ht::MonomialHashtable)
    @inbounds e1 = ht.exponents[h1]
    @inbounds e2 = ht.exponents[h2]
    is_gcd_const(e1, e2)
end

# computes lcm of he1 and he2 as exponent vectors from ht1
# and inserts it in ht2
function get_lcm(he1::MonomIdx, he2::MonomIdx,
    ht1::MonomialHashtable{M}, ht2::MonomialHashtable{M}) where {M}

    @inbounds e1 = ht1.exponents[he1]
    @inbounds e2 = ht1.exponents[he2]
    @inbounds etmp = ht1.exponents[1]

    etmp = monom_lcm!(etmp, e1, e2)

    insert_in_hash_table!(ht2, etmp)
end

#------------------------------------------------------------------------------

# compare pairwise divisibility of lcms from a[first:last] with lcm
function check_monomial_division_in_update(
    a::Vector{MonomIdx}, first::Int, last::Int,
    lcm::MonomIdx, ht::MonomialHashtable{M}) where {M}

    # pairs are sorted, we only need to check entries above starting point

    @inbounds divmask = ht.hashdata[lcm].divmask
    @inbounds lcmexp = ht.exponents[lcm]

    j = first
    @inbounds while j <= last
        # bad lcm
        if iszero(a[j])
            j += 1
            continue
        end
        # fast division check
        if !is_divmask_divisible(ht.hashdata[a[j]].divmask, divmask)
            j += 1
            continue
        end
        ea = ht.exponents[a[j]]
        if !is_monom_divisible(ea, lcmexp)
            j += 1
            continue
        end
        # mark as redundant
        a[j] = 0
    end

    nothing
end

#------------------------------------------------------------------------------

# add monomials from `poly` multiplied by exponent vector `etmp`
# with hash `htmp` to hashtable `symbol_ht`,
# and substitute hashes in row
function insert_multiplied_poly_in_hash_table!(
        row::Vector{MonomIdx},
        htmp::MonomHash,
        etmp::M,
        poly::Vector{MonomIdx},
        ht::MonomialHashtable{M},
        symbol_ht::MonomialHashtable{M}) where {M}

    # length of poly to add
    len = length(poly)

    mod = MonomHash(symbol_ht.size - 1)

    bexps = ht.exponents
    bdata = ht.hashdata

    sexps = symbol_ht.exponents
    sdata = symbol_ht.hashdata

    l = 1 # hardcoding 1 does not seem nice =(
    @label Letsgo
    @inbounds while l <= len
        # we iterate over all monoms of the given poly,
        # multiplying them by htmp/etmp,
        # and inserting into symbolic hashtable

        # hash is linear, so that
        # hash(e1 + e2) = hash(e1) + hash(e2)
        # We also assume that the hashing vector is shared same
        # between all created hashtables
        h = htmp + bdata[poly[l]].hash

        e = bexps[poly[l]]

        lastidx = symbol_ht.load + 1
        enew = sexps[1]
        enew = monom_product!(enew, etmp, e)

        # insert into hashtable
        k = h

        i = MonomHash(1)
        ssize = MonomHash(symbol_ht.size)
        @inbounds while i <= ssize
            k = nexthashindex(h, i, mod)
            vidx = symbol_ht.hashtable[k]

            # if index is free
            iszero(vidx) && break

            if ishashcollision(symbol_ht, vidx, enew, h)
                i += MonomHash(1)
                continue
            end
            
            # hit
            row[l] = vidx
            l += 1

            @goto Letsgo
        end
        # miss

        # add multiplied exponent to hash table        
        sexps[lastidx] = copy(enew)
        symbol_ht.hashtable[k] = lastidx

        divmask = monom_divmask(enew, DivisionMask, 
                    symbol_ht.ndivvars, symbol_ht.divmap, symbol_ht.ndivbits)
        sdata[lastidx] = Hashvalue(0, h, divmask, totaldeg(enew))

        row[l] = lastidx
        l += 1
        symbol_ht.load += 1
    end

    row
end

function multiplied_poly_to_matrix_row!(
    symbolic_ht::MonomialHashtable, basis_ht::MonomialHashtable{M},
    htmp::MonomHash, etmp::M, poly::Vector{MonomIdx}) where {M}

    row = similar(poly)
    check_enlarge_hashtable!(symbolic_ht, length(poly))

    insert_multiplied_poly_in_hash_table!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
end

#------------------------------------------------------------------------------

function insert_in_basis_hash_table_pivots(
    row::Vector{ColumnIdx},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M},
    col2hash::Vector{MonomIdx}) where {M}

    check_enlarge_hashtable!(ht, length(row))

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.exponents

    mod = MonomHash(ht.size - 1)
    bdata = ht.hashdata
    bexps = ht.exponents
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        hidx = col2hash[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = MonomHash(1)
        @inbounds while i <= ht.size
            k = nexthashindex(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if ishashcollision(ht, hm, e, h)
                i += MonomHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(sdata[hidx].idx, h, sdata[hidx].divmask, sdata[hidx].deg)

        ht.load += 1
    end

    nothing
end
