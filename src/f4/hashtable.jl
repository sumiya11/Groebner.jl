# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

### 
# Monomial hashtable.

# Some monomial implementations are mutable and some are not. To
# maintain generic code that will work for both, we usually write
#   m3 = monom_product!(m3, m1, m2). 
# If a mutable implementation is used, m3 would be overwritten inside
# monom_product!. Otherwise, monom_product! would return a new immutable object
# that is then assigned to m3. That allows us to write more or less
# independently of the monomial implementation.

# Hash of a monomial in the hashtable
const MonomHash = UInt32

# The idenfifier of a monomial. This idenfifier is guaranteed to be unique
# within a particular hashtable. This allows one to use this idenfifier when
# working with monomials
const MonomId = Int32

# Division mask of a monomial
const DivisionMask = UInt32

# Hashtable implements open addressing with linear scan.
mutable struct MonomialHashtable{M <: Monom, Ord <: AbstractMonomialOrdering}
    #= Data =#
    monoms::Vector{M}
    hashtable::Vector{MonomId}
    labels::Vector{Int32}
    hashvals::Vector{MonomHash}
    divmasks::Vector{DivisionMask}

    hash_vector::Vector{MonomHash}

    #= Ring information =#
    nvars::Int
    ord::Ord

    #= Monom divisibility =#
    use_divmask::Bool
    compress_divmask::Bool
    divmap::Vector{UInt32}
    ndivvars::Int
    ndivbits::Int

    # Hashtable sizes 
    size::Int
    load::Int
    offset::Int

    # If the hashtable is frozen, modification is not allowed
    frozen::Bool
end

###
# Initialization

hashtable_resize_threshold() = 0.4
hashtable_needs_resize(size, load, added) = (load + added) / size > hashtable_resize_threshold()

function hashtable_initialize(
    ring::PolyRing{Ord},
    rng::AbstractRNG,
    MonomT::T
) where {Ord <: AbstractMonomialOrdering, T}
    initial_size = 2^10

    hashtable = zeros(MonomId, initial_size)

    # The table is filled to no more than 50% at any time so there is no need to
    # allocate the full storage
    half_size = div(initial_size, 2)
    @assert half_size > 0
    monoms = Vector{MonomT}(undef, half_size)
    labels = Vector{Int32}(undef, half_size)
    hashvals = Vector{MonomHash}(undef, half_size)
    divmasks = Vector{DivisionMask}(undef, half_size)

    nvars = ring.nvars
    ord = ring.ord

    hash_vector = monom_construct_hash_vector(rng, MonomT, nvars)

    load = 1
    @invariant initial_size > 1
    size = initial_size

    # Exponents array starts from index 2,
    # We store buffer array at index 1.
    offset = 2

    # Initialize fast divisibility parameters.
    use_divmask = nvars <= 32 || (MonomT <: ExponentVector)
    dmbits = 8 * sizeof(DivisionMask)
    compress_divmask = nvars > dmbits
    ndivbits = div(dmbits, nvars)
    # division mask stores at least 1 bit
    # per each of first ndivvars variables
    ndivbits == 0 && (ndivbits += 1)
    # count only first ndivvars variables for divisibility checks
    ndivvars = min(nvars, dmbits)
    divmap = zeros(DivisionMask, ndivvars * ndivbits)

    # first stored exponent used as buffer lately
    monoms[1] = monom_construct_const(MonomT, nvars)

    MonomialHashtable(
        monoms,
        hashtable,
        labels,
        hashvals,
        divmasks,
        hash_vector,
        nvars,
        ord,
        use_divmask,
        compress_divmask,
        divmap,
        ndivvars,
        ndivbits,
        size,
        load,
        offset,
        false
    )
end

function hashtable_initialize_secondary(ht::MonomialHashtable{M}) where {M <: Monom}
    # 2^6 seems to be the best out of 2^5, 2^6, 2^7
    initial_size = 2^6
    @invariant initial_size > 1

    hashtable = zeros(MonomId, initial_size)

    half_size = div(initial_size, 2)
    monoms = Vector{M}(undef, half_size)
    labels = Vector{Int32}(undef, half_size)
    hashvals = Vector{MonomHash}(undef, half_size)
    divmasks = Vector{DivisionMask}(undef, half_size)

    # preserve ring info
    nvars = ht.nvars
    ord = ht.ord

    # preserve division info
    divmap = ht.divmap
    ndivbits = ht.ndivbits
    ndivvars = ht.ndivvars

    # preserve hash_vector
    hash_vector = ht.hash_vector

    load = 1
    size = initial_size
    offset = 2

    monoms[1] = monom_construct_const(M, nvars)

    MonomialHashtable(
        monoms,
        hashtable,
        labels,
        hashvals,
        divmasks,
        hash_vector,
        nvars,
        ord,
        ht.use_divmask,
        ht.compress_divmask,
        divmap,
        ndivvars,
        ndivbits,
        size,
        load,
        offset,
        false
    )
end

function hashtable_reinitialize!(ht::MonomialHashtable{M}) where {M}
    @invariant !ht.frozen

    initial_size = 2^6
    @invariant initial_size > 1

    # Reinitialize counters
    ht.load = 1
    ht.offset = 2
    ht.size = initial_size

    resize!(ht.hashtable, ht.size)

    half_size = div(ht.size, 2)
    resize!(ht.monoms, half_size)
    resize!(ht.labels, half_size)
    resize!(ht.hashvals, half_size)
    resize!(ht.divmasks, half_size)

    hashtable = ht.hashtable
    @inbounds for i in 1:(ht.size)
        hashtable[i] = zero(MonomId)
    end

    ht.monoms[1] = monom_construct_const(M, ht.nvars)

    nothing
end

function hashtable_resize_if_needed!(ht::MonomialHashtable, added::Int)
    newsize = ht.size
    while hashtable_needs_resize(newsize, ht.load, added)
        newsize *= 2
    end
    newsize == ht.size && return nothing

    @invariant !ht.frozen
    @invariant ispow2(newsize)

    ht.size = newsize
    resize!(ht.hashtable, ht.size)

    half_size = div(ht.size, 2)
    @assert half_size >= ht.load
    resize!(ht.monoms, half_size)
    resize!(ht.labels, half_size)
    resize!(ht.hashvals, half_size)
    resize!(ht.divmasks, half_size)

    @inbounds for i in 1:(ht.size)
        ht.hashtable[i] = zero(MonomId)
    end

    mod = MonomHash(ht.size - 1)

    @inbounds for i in (ht.offset):(ht.load)
        # hash for this elem is already computed
        he = ht.hashvals[i]
        hidx = he
        for j in MonomHash(0):MonomHash(ht.size)
            hidx = hashtable_next_lookup_index(he, j, mod)
            !iszero(ht.hashtable[hidx]) && continue
            ht.hashtable[hidx] = i
            break
        end
    end
    nothing
end

###
# Insertion of monomials

# Must be within 1 <= ... <= mod+1
function hashtable_next_lookup_index(h::MonomHash, j::MonomHash, mod::MonomHash)
    ((h + j) & mod) + MonomHash(1)
end

function hashtable_is_hash_collision(ht::MonomialHashtable, vidx, e, he)
    @inbounds ht.hashvals[vidx] != he && return true
    @inbounds !monom_is_equal(ht.monoms[vidx], e) && return true
    false
end

function hashtable_insert!(ht::MonomialHashtable{M}, e::M) where {M <: Monom}
    # Optimizing for the case when the monomial is already in the table.
    # All of the functions in the main code path are inlined.
    @invariant ispow2(ht.size) && ht.size > 1

    he = monom_hash(e, ht.hash_vector)

    hsize = ht.size
    mod = (hsize - 1) % MonomHash
    hidx = hashtable_next_lookup_index(he, 0 % MonomHash, mod)
    @inbounds vidx = ht.hashtable[hidx]

    hit = !iszero(vidx)
    @inbounds if hit && !hashtable_is_hash_collision(ht, vidx, e, he)
        # Hit!
        return vidx
    end

    # Miss or collision
    i = 1 % MonomHash
    mhhsize = hsize % MonomHash
    @inbounds while hit && i < mhhsize
        hidx = hashtable_next_lookup_index(he, i, mod)
        vidx = ht.hashtable[hidx]

        iszero(vidx) && break

        if hashtable_is_hash_collision(ht, vidx, e, he)
            i += (1 % MonomHash)
            continue
        end

        # already present in hashtable
        return vidx
    end

    @invariant !ht.frozen

    # add monomial to hashtable
    vidx = (ht.load + 1) % MonomId
    @invariant vidx <= length(ht.monoms)
    @inbounds ht.hashtable[hidx] = vidx
    @inbounds ht.monoms[vidx] = monom_copy(e)
    divmask = monom_create_divmask(
        e,
        DivisionMask,
        ht.ndivvars,
        ht.divmap,
        ht.ndivbits,
        ht.compress_divmask
    )
    @inbounds ht.labels[vidx] = NON_PIVOT_COLUMN
    @inbounds ht.hashvals[vidx] = he
    @inbounds ht.divmasks[vidx] = divmask

    ht.load += 1

    return vidx
end

###
# Division masks

function divmask_is_probably_divisible(d1::DivisionMask, d2::DivisionMask)
    iszero(~d1 & d2)
end

function hashtable_fill_divmasks!(ht::MonomialHashtable)
    @invariant !ht.frozen

    ndivvars = ht.ndivvars

    min_exp = Vector{UInt64}(undef, ndivvars)
    max_exp = Vector{UInt64}(undef, ndivvars)

    e = Vector{UInt64}(undef, ht.nvars)
    monom_to_vector!(e, ht.monoms[ht.offset])

    @inbounds for i in 1:ndivvars
        min_exp[i] = e[i]
        max_exp[i] = e[i]
    end

    @inbounds for i in (ht.offset):(ht.load)
        monom_to_vector!(e, ht.monoms[i])
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

    if !ht.compress_divmask
        # Available bits >= variables
        ctr = 1
        steps = UInt32(0)
        @inbounds for i in 1:ndivvars
            steps = div(max_exp[i] - min_exp[i], UInt32(ht.ndivbits))
            (iszero(steps)) && (steps += UInt32(1))
            for j in 1:(ht.ndivbits)
                ht.divmap[ctr] = steps
                steps += UInt32(1)
                ctr += 1
            end
        end
    else
        # Available bits < variables.
        # Pack variables tighlty.
        vars_per_bit = div(ht.nvars, 8 * sizeof(DivisionMask))
        vars_covered = vars_per_bit * 8 * sizeof(DivisionMask)
        if vars_covered != ht.nvars
            @invariant vars_covered < ht.nvars
            vars_per_bit += 1
        end
        @invariant vars_per_bit > 1
        @invariant ht.ndivbits == 1
        @invariant ndivvars == length(ht.divmap)
        ctr = 1
        @inbounds for i in 1:ndivvars
            if ht.nvars - ctr + 1 <= (ndivvars - i + 1) * (vars_per_bit - 1)
                vars_per_bit -= 1
            end
            ht.divmap[i] = vars_per_bit
            ctr += vars_per_bit
        end
        @invariant sum(ht.divmap) == ht.nvars
    end

    @inbounds for vidx in (ht.offset):(ht.load)
        divmask = monom_create_divmask(
            ht.monoms[vidx],
            DivisionMask,
            ht.ndivvars,
            ht.divmap,
            ht.ndivbits,
            ht.compress_divmask
        )
        ht.labels[vidx] = NON_PIVOT_COLUMN
        ht.divmasks[vidx] = divmask
    end

    nothing
end

###
# Monomial arithmetic

function hashtable_monom_is_divisible(h1::MonomId, h2::MonomId, ht::MonomialHashtable)
    @inbounds if ht.use_divmask
        if !divmask_is_probably_divisible(ht.divmasks[h1], ht.divmasks[h2])
            return false
        end
    end
    @inbounds e1 = ht.monoms[h1]
    @inbounds e2 = ht.monoms[h2]
    monom_is_divisible(e1, e2)
end

function hashtable_get_lcm!(
    he1::MonomId,
    he2::MonomId,
    ht1::MonomialHashtable{M},
    ht2::MonomialHashtable{M}
) where {M}
    @inbounds e1 = ht1.monoms[he1]
    @inbounds e2 = ht1.monoms[he2]
    @inbounds etmp = ht1.monoms[1]

    etmp = monom_lcm!(etmp, e1, e2)

    hashtable_insert!(ht2, etmp)
end

###

# Inserts a multiple of the polynomial into symbolic hashtable.
# Writes the resulting monomial identifiers to the given row.
function hashtable_insert_polynomial_multiple!(
    row::Vector{MonomId},
    mult_hash::MonomHash,
    mult::M,
    poly::Vector{MonomId},
    ht::MonomialHashtable{M},
    symbol_ht::MonomialHashtable{M}
) where {M <: Monom}
    @invariant ispow2(ht.size) && ht.size > 1
    @invariant ispow2(symbol_ht.size) && symbol_ht.size > 1
    @invariant length(row) == length(poly)

    iszero(length(poly)) && return row

    ssize = symbol_ht.size % MonomHash
    mod = (symbol_ht.size - 1) % MonomHash
    @inbounds buf = symbol_ht.monoms[1]

    # Iterate over monomials of the given polynomial, multiply them by a
    # monomial multiple, and insert them into symbolic hashtable. 
    # We use the fact that the hash function is linear.
    #
    # It is often the case that the multiple of the leading monomial is already
    # in the symbolic hashtable. In this case, skip it.
    @inbounds for j in 1:length(poly)
        oldmonom = ht.monoms[poly[j]]
        newmonom = monom_product!(buf, mult, oldmonom)

        oldhash = ht.hashvals[poly[j]]
        newhash = mult_hash + oldhash

        hidx = hashtable_next_lookup_index(newhash, 0 % MonomHash, mod)
        vidx = symbol_ht.hashtable[hidx]

        hit = !iszero(vidx)
        if hit && !hashtable_is_hash_collision(symbol_ht, vidx, newmonom, newhash)
            # Hit, go to next monomial.
            row[j] = vidx
            continue
        end

        # Miss or collision.
        i = 1 % MonomHash
        while hit && i <= ssize
            hidx = hashtable_next_lookup_index(newhash, i, mod)
            vidx = symbol_ht.hashtable[hidx]

            # Found a free bucket.
            hit = !iszero(vidx)
            iszero(vidx) && break

            if hashtable_is_hash_collision(symbol_ht, vidx, newmonom, newhash)
                i += (1 % MonomHash)
                continue
            end

            # Hit, go to next monomial.
            row[j] = vidx
            break
        end
        hit && continue
        # Miss! Add monomial multiple to the hash table.

        @invariant !symbol_ht.frozen

        vidx = (symbol_ht.load + 1) % MonomId
        @invariant vidx <= length(symbol_ht.monoms)
        symbol_ht.monoms[vidx] = monom_copy(newmonom)
        symbol_ht.hashtable[hidx] = vidx
        divmask = monom_create_divmask(
            newmonom,
            DivisionMask,
            symbol_ht.ndivvars,
            symbol_ht.divmap,
            symbol_ht.ndivbits,
            symbol_ht.compress_divmask
        )
        @inbounds symbol_ht.labels[vidx] = NON_PIVOT_COLUMN
        @inbounds symbol_ht.hashvals[vidx] = newhash
        @inbounds symbol_ht.divmasks[vidx] = divmask

        row[j] = vidx
        symbol_ht.load += 1
    end

    row
end
