
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
####### Hashtable ######

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

function hashnextindex(h::UInt32, j::UInt32, mod::UInt32)
     (h + j) & mod + UInt32(1)
end

#------------------------------------------------------------------------------

# initialize and set fields for basis hashtable
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

function copy_hashtable(ht::MonomialHashtable)
    exps = Vector{Vector{UInt16}}(undef, ht.size)
    table = Vector{Int}(undef, ht.size)
    data = Vector{Hashvalue}(undef, ht.size)
    exps[1] = zeros(UInt16, ht.explen)

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

    # preserve hasher
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
####### Basis #######

# s-pair, a pair of polynomials
struct SPair
    # first generator as index from the basis array
    poly1::Int
    # second generator -//-
    poly2::Int
    # position of lcm(poly1, poly2) in hashtable
    lcm::Int
    # total degree of lcm
    deg::UInt
end

mutable struct Pairset
    pairs::Vector{SPair}
    # number of filled pairs,
    # Initially zero
    load::Int
end

function initialize_pairset(; initial_size=2^6) # TODO: why 64?
    pairs = Vector{SPair}(undef, initial_size)
    return Pairset(pairs, 0)
end

function Base.isempty(ps::Pairset)
    return ps.load == 0
end

function check_enlarge_pairset!(ps::Pairset, added::Int)
    sz = length(ps.pairs)
    # TODO: shrink pairset ?
    if ps.load + added >= sz
        newsz = max(2 * sz, ps.load + added)
        resize!(ps.pairs, newsz)
    end
end

#------------------------------------------------------------------------------

mutable struct Basis{T}
    # vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with it's position in hashtable
    gens::Vector{Vector{Int}}
    # polynomial coefficients
    coeffs::Vector{Vector{T}}

    #= Keeping track of sizes   =#
    #=  ndone <= ntotal <= size =#
    # total allocated size,
    # size == length(gens) is always true
    size::Int
    # number of processed polynomials in `gens`
    # Iitially a zero
    ndone::Int
    # total number of polys filled in `gens`
    # (these will be handled during next update! call)
    # Iitially a zero
    ntotal::Int

    #= Keeping track of redundancy =#
    #= invariant =#
    #= length(lead) == length(nonred) == count(isred) == nlead =#
    # if element of the basis
    # is redundant
    isred::Vector{Int8}
    # positions of non-redundant elements in the basis
    nonred::Vector{Int}
    # division masks of leading monomials of
    # non redundant basis elements
    lead::Vector{UInt32}
    # number of filled elements in lead
    nlead::Int

    # characteristic of ground field
    ch::UInt64
end

function initialize_basis(ring::PolyRing, ngens::Int, ::Type{T}) where {T<:Coeff}
    #=
        always true
        length(gens) == length(coeffs) == length(isred) == size
    =#

    sz     = ngens * 2 # hmmm
    ndone  = 0
    ntotal = 0
    nlead  = 0

    gens   = Vector{Vector{Int}}(undef, sz)
    coeffs = Vector{Vector{T}}(undef, sz)
    isred  = zeros(Int8, sz)
    nonred = Vector{Int}(undef, sz)
    lead   = Vector{UInt32}(undef, sz)

    ch = ring.ch

    Basis(gens, coeffs, sz, ndone, ntotal, isred, nonred, lead, nlead, ch)
end

function initialize_basis(ring::PolyRing, hashedexps, coeffs::Vector{Vector{T}}) where {T<:Coeff}
    sz     = length(hashedexps) # hmmm
    ndone  = 0
    ntotal = 0
    nlead  = 0

    isred  = zeros(Int8, sz)
    nonred = Vector{Int}(undef, sz)
    lead   = Vector{UInt32}(undef, sz)

    ch = ring.ch

    Basis(hashedexps, coeffs, sz, ndone, ntotal, isred, nonred, lead, nlead, ch)
end

#------------------------------------------------------------------------------

function copy_basis_thorough(basis::Basis{T}) where {T}
    #  That cost a day of debugging ////
    gens   = Vector{Vector{Int}}(undef, basis.size)
    coeffs = Vector{Vector{T}}(undef, basis.size)
    @inbounds for i in 1:basis.ntotal
        gens[i] = Vector{Int}(undef, length(basis.gens[i]))
        coeffs[i] = Vector{T}(undef, length(basis.coeffs[i]))

        @inbounds for j in 1:length(basis.gens[i])
            gens[i][j] = basis.gens[i][j]
            coeffs[i][j] = basis.coeffs[i][j]
        end
    end
    isred  = copy(basis.isred)
    nonred = copy(basis.nonred)
    lead = copy(basis.lead)
    Basis(gens, coeffs, basis.size, basis.ndone,
            basis.ntotal, isred, nonred, lead,
            basis.nlead, basis.ch)
end

#------------------------------------------------------------------------------

function check_enlarge_basis!(basis::Basis{T}, added::Int) where {T}
    if basis.ndone + added >= basis.size
        basis.size = max(basis.size * 2, basis.ndone + added)
        resize!(basis.gens, basis.size)
        resize!(basis.coeffs, basis.size)
        resize!(basis.isred, basis.size)
        basis.isred[basis.ndone+1:end] .= 0
        resize!(basis.nonred, basis.size)
        resize!(basis.lead, basis.size)
    end
end

#------------------------------------------------------------------------------

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(basis::Basis{CoeffFF})
    cfs = basis.coeffs
    @inbounds for i in 1:basis.ntotal
        # mul = inv(cfs[i][1])
        # hack for now, TODODO
        if !isassigned(cfs, i)
            continue
        end
        ch = basis.ch
        mul = invmod(cfs[i][1], ch) % ch
        @inbounds for j in 2:length(cfs[i])
            # cfs[i][j] *= mul
            cfs[i][j] = cfs[i][j]*mul % ch
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(basis::Basis{CoeffQQ})
    cfs = basis.coeffs
    @inbounds for i in 1:basis.ntotal
        # mul = inv(cfs[i][1])
        # hack for now, TODODO
        if !isassigned(cfs, i)
            continue
        end
        mul = inv(cfs[i][1])
        @inbounds for j in 2:length(cfs[i])
            cfs[i][j] *= mul
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

#------------------------------------------------------------------------------
####### Matrix #######

mutable struct MacaulayMatrix{T<:Coeff}
    #=
        Matrix of the following structure

        | A  B |
        | C  D |

        A contains known pivots of reducing rows,
        and CD are rows to be reduced by AB

    =#

    # rows from upper, AB part of the matrix,
    # stored as vectors of corresponding exponents (already hashed)
    uprows::Vector{Vector{Int}}
    # rows from lower, CD part of the matrix,
    # stored as vectors of corresponding exponents (already hashed)
    lowrows::Vector{Vector{Int}}

    # maps column idx {1 ... ncols} to monomial hash {2 ... ht.load}
    # in some (?) hashtable
    col2hash::Vector{Int}

    # row coefficients
    # (some of the rows are stored in the basis,
    #  and some are stored here)
    coeffs::Vector{Vector{T}}

    #= sizes info =#
    # total number of allocated rows
    size::Int
    # number of pivots,
    # ie new basis elements discovered after matrix reduction
    npivots::Int
    # number of filled rows, nrows <= size
    nrows::Int
    # number of columns
    ncols::Int

    # number of upper rows (in AB section)
    nup::Int
    # number of lower rows (in CD section)
    nlow::Int
    # number of left cols  (in AC section)
    nleft::Int
    # number of right cols (in BD section)
    nright::Int

    # maps column idx {1 ... ncols} to index of coefficient array
    # First nleft indices point to coefficients from AB part of the matrix
    # (these coefficients are owned by the basis struct)
    # Last nright indices point to coefficients from CD part of the matrix
    # (these are owned by this object and stored in coeffs array)
    # Essentially each coefficient array from the first part represents a reducer row,
    # and each array from the second part stands for a reduced *nonzero* row
    # should be row to coef
    up2coef::Vector{Int}
    low2coef::Vector{Int}
end

function initialize_matrix(ring::PolyRing, ::Type{T}) where {T<:Coeff}
    uprows   = Vector{Vector{Int}}(undef, 0)
    lowrows  = Vector{Vector{Int}}(undef, 0)
    col2hash = Vector{Int}(undef, 0)
    coeffs   = Vector{Vector{T}}(undef, 0)

    size    = 0
    npivots = 0
    nrows   = 0
    ncols   = 0

    nup    = 0
    nlow   = 0
    nleft  = 0
    nright = 0

    up2coef   = Vector{Int}(undef, 0)
    low2coef   = Vector{Int}(undef, 0)

    MacaulayMatrix(uprows, lowrows, col2hash, coeffs,
            size, npivots, nrows, ncols,
            nup, nlow, nleft, nright,
            up2coef, low2coef)
end

# TODO: change semantic
function reinitialize_matrix!(matrix::MacaulayMatrix{T}, npairs::Int) where {T}
    resize!(matrix.uprows, npairs*2)
    resize!(matrix.lowrows, npairs*2)
    resize!(matrix.up2coef, npairs*2)
    resize!(matrix.low2coef, npairs*2)
    matrix.size = 2*npairs
    matrix.ncols = 0
    matrix.nleft = 0
    matrix.nright = 0
    matrix.nup   = 0
    matrix.nlow  = 0
    matrix
end
