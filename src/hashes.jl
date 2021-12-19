
# todo: implementat custom hashtable with linear probing

# Hashtable for monomials occurring during computations
# Allows us to operate monomials as int hashes at the upper level
#
# Stores exponents in raw format as in AA
struct MonomialHashtable{Tv}
    # maps hash, an int, to corresponding exponent vector,
    # HashedMonomial objects wrap hash
    expmap::Dict{HashedMonomial, Vector{UInt}}
    # base ring
    R::MPolyRing{Tv}
    # values to hash with
    hashvals::Vector{UInt}

    # length of exponent vector (+1 if degrevlex)
    explen::UInt
end

function MonomialHashtable(R::MPolyRing{Tv}) where {Tv}
    expmap = Dict{HashedMonomial, Vector{UInt}}()
    Base.rehash!(expmap, 2^16)  # hmmm

    explen = nvars(R)
    if ordering(R) == :degrevlex
        explen += 1
    end

    hashvals = rand(UInt, nvars(R))

    MonomialHashtable{Tv}(expmap, R, hashvals, UInt(explen))
end

#------------------------------------------------------------------------------

# adds monomial represented with exponent vector to hash table
# and returns a new hash
function add_monomial!(HT::MonomialHashtable, exp::Vector{UInt})
    # calculate exponent hash
    hexp = UInt(0)
    for i in 1:HT.explen
        hexp += HT.hashvals[i] * exp[i]
    end

    hm = HashedMonomial(hexp)
    (haskey(HT.expmap, hm)) && (return hm)

    HT.expmap[hm] = exp
    hm
end

#=
function add_product!(HT::MonomialHashtable, h1, h2)
    (haskey(HT.prods, (h1, h2))) && (return HT.prods[(h1, h2)])

    # make hash additive
    prode = HT.expmap[h1] .+ HT.expmap[h2]
    prodh = add_monomial!(HT, prode)

    HT.prods[(h1, h2)] = prodh

    prodh
end
=#

# returns exponent vector for the given hash
function Base.getindex(HT::MonomialHashtable, h::HashedMonomial)
    getindex(HT.expmap, h)
end
