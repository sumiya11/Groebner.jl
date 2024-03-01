
# Hash of a monomial in the hashtable
const MonomHash = UInt32

# Hashvalue1 of a single monomial
struct Hashvalue1
    # index of the monomial in the F4 matrix (defaults to NON_PIVOT_COLUMN, or 0),
    idx::Int32
    # hash of the monomial,
    hash::MonomHash
    # total degree of the monomial
    deg::MonomHash
end

# Hashtable implements open addressing with linear scan.
mutable struct MonomialHashtable1{M <: Monom}
    #= Data =#
    monoms::Vector{M}
    # Maps monomial id to its position in the `monoms` array
    hashtable::Vector{MonomId}
    # Stores hashes, division masks, and other valuable info for each hashtable
    # enrty
    hashdata::Vector{Hashvalue1}
    # Hash vector. Hash of a monomial is a dot product of the `hasher` vector
    # and the monomial exponent vector
    hasher::Vector{MonomHash}

    #= Ring information =#
    # number of variables
    nvars::Int

    # Hashtable size 
    # (always a power of two)
    # (always greater than 1)
    size::Int
    # Elements currently added
    load::Int
    offset::Int

    # If the hashtable is frozen, any operation that tries to modify it will
    # result in an error.
    frozen::Bool
end

###
# Initialization and resizing

# Resize hashtable if load factor exceeds hashtable_resize_threshold. Load factor of a
# hashtable must be smaller than hashtable_resize_threshold at any point of its
# lifetime
hashtable_resize_threshold() = 0.4
hashtable_needs_resize(size, load, added) =
    (load + added) / size > hashtable_resize_threshold()

function hashtable_initialize1(
    nvars,
    rng::AbstractRNG,
    MonomT::T,
    initial_size::Int
) where {T}
    exponents = Vector{MonomT}(undef, initial_size)
    hashdata = Vector{Hashvalue1}(undef, initial_size)
    hashtable = zeros(MonomId, initial_size)

    # initialize hashing vector
    hasher = [rand(MonomHash) for i in 1:nvars]

    # exponents[1:load] covers all stored exponents
    # , also exponents[1] is [0, 0, ..., 0] by default
    load = 1
    @assert initial_size > 1
    size = initial_size

    # exponents array starts from index offset,
    # We store buffer array at index 1
    offset = 2

    # first stored exponent used as buffer lately
    exponents[1] = zeros(UInt8, nvars)

    MonomialHashtable1(
        exponents,
        hashtable,
        hashdata,
        hasher,
        nvars,
        size,
        load,
        offset,
        false
    )
end

function hashtable_resize_if_needed!(ht::MonomialHashtable1, added::Int)
    newsize = ht.size
    while hashtable_needs_resize(newsize, ht.load, added)
        newsize *= 2
    end
    newsize == ht.size && return nothing

    ht.size = newsize

    resize!(ht.hashdata, ht.size)
    resize!(ht.monoms, ht.size)
    resize!(ht.hashtable, ht.size)
    @inbounds for i in 1:(ht.size)
        ht.hashtable[i] = zero(MonomId)
    end

    mod = MonomHash(ht.size - 1)

    @inbounds for i in (ht.offset):(ht.load)
        # hash for this elem is already computed
        he = ht.hashdata[i].hash
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

# Returns the next look-up position in the table.
# Must be within 1 <= ... <= mod+1
function hashtable_next_lookup_index(h::MonomHash, j::MonomHash, mod::MonomHash)
    ((h + j) & mod) + MonomHash(1)
end

# if hash collision happened
function hashtable_is_hash_collision(ht::MonomialHashtable1, vidx, e, he)
    # if not free and not same hash
    @inbounds if ht.hashdata[vidx].hash != he
        return true
    end
    # if not free and not same monomial
    @inbounds if !(ht.monoms[vidx] == e)
        return true
    end
    false
end

function monom_hash(x::Vector{T}, b::Vector{MH}) where {T, MH}
    h = zero(MH)
    @inbounds for i in eachindex(x, b)
        h = h + MH(x[i]) * b[i]
    end
    mod(h, MonomHash)
end

function hashtable_insert!(ht::MonomialHashtable1{M}, e::M) where {M <: Monom}
    # NOTE: trying to optimize for the case when the monomial is already in the
    # table.
    # NOTE: all of the functions called here are inlined. The only potential
    # exception is monom_create_divmask

    # generate hash
    he = monom_hash(e, ht.hasher)

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

    # add monomial to hashtable
    vidx = (ht.load + 1) % MonomId
    @inbounds ht.hashtable[hidx] = vidx
    @inbounds ht.monoms[vidx] = Base.copy(e)
    @inbounds ht.hashdata[vidx] = Hashvalue1(0, he, e[1])

    ht.load += 1

    return vidx
end
