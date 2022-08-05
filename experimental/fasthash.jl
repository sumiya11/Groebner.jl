
#=
    Hashtable:
    
    a Monomial is an Integer

    Integer is an index to hashtable, which maps integer -> exponent vector

    hashtable has two parts:
    
=#

#=

degree 1 = (degree 1)^0
    [0,.. 0]

degree 1 = degree 1
    [1,0..0], [0,1,0...0], ..[0,...,0,1]

degree 2 = degree 1 × degree 1
    ...

degree 3 = degree 2 × degree 1
    ...

degree 4 = degree 3 × degree 1
    ...

...

=#


import Random

function nallterms(n, d)
    d == 0 && return BigInt(1)
    nallterms(n, d - 1) + BigInt(n)^d
end

const ExponentVector = Vector{Int}

#=
    stores hash of a monomial,
    corresponding divmask to speed up divisibility checks,
    index for position matrix (defaults to zero),
    and the todal degree
=#
mutable struct Hashvalue
    hash::Int
    deg::Int
end

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
    hasher::Vector{Int}

    #= Ring information =#
    # number of variables
    nvars::Int

    size::Int
    # elements added
    load::Int
    # 
    offset::Int
end


# Returns the next look-up position in the table 
# (that is, implementing open addressing with quadratic probing) 
function hashnextindex(h::Int, j::Int, mod::Int)
    (h + j) & mod + 1
end

function ithvector(nvars, deg, i)
    # nvars = 3
    # deg = 2
    #=
    0 [0, 0, 0]  + 0
    1 [0, 0, 1]
    2 [0, 0, 2]
    3 [0, 1, 0]  i - (deg+1)^1
    4 [0, 1, 1]
    5 [0, 1, 2]
    6 [0, 2, 0]
    7 [0, 2, 1]
    8 [0, 2, 2]
    9 [1, 0, 0]  i - (deg + 1)^2
    10 [1, 0, 1]  
    11 [1, 0, 2]
    12 [1, 1, 0]
    13 [1, 1, 1]
    14 [1, 1, 2]
    15 [1, 2, 0]
    ...
    =#
    s = ExponentVector(undef, nvars+1)
    p = nvars-1
    for j in 2:nvars+1
        s[j] = div(i, (deg+1)^p)
        i -= div(i, (deg+1)^p)*(deg+1)^p
        p -= 1
    end
    s[1] = sum((s[j] for j in 2:nvars+1))
    s
end

function monomprod(ht, i1, i2)
    if i1 < ht.offset
        m1 = ht.exponents[i1]
    end
    if i2 < ht.offset
        m2 = ht.exponents[i2]
    end
    
    # i1 * i2 = 1

    # if 
end

# initialize and set fields for basis hashtable
function initialize_hash_table(
    nvars::Int, deg::Int;
    cache=-1,
    rng::Random.AbstractRNG=Random.MersenneTwister(),
    initial_size::Int=2^3)

    exponents = Vector{ExponentVector}(undef, initial_size)
    hashdata  = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(Int, initial_size)

    nvars = nvars

    # initialize hashing vector
    hasher = [Int(i) for i in 1:nvars+1]
    
    # exponents[1:load] cover all stored exponents
    # , also exponents[1] is zeroed by default
    load = 1
    size = initial_size

    offset = (deg+1)^nvars
    if cache != -1
        offset = cache+1
    end

    for i in 1:offset
        exponents[i] = ithvector(nvars, deg, i - 1) 
        hashtable[i] = i
    end

    MonomialHashtable(
        exponents, hashtable, hashdata, hasher,
        nvars, size, load, offset)
end

function vars(ht)
    1, 2, 3
end

ht = initialize_hash_table(2, 2, initial_size=2^4)

one, x1, x2 = vars(ht)

monomprod(ht, one, x1)