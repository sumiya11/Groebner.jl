
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
function add_monomial!(ht::MonomialHashtable, exp::Vector{UInt})
    # calculate exponent hash
    #=
    hexp = UInt(0)
    for i in 1:ht.explen
        hexp += ht.hashvals[i] * exp[i]
    end
    =#

    hexp = hash(exp)
    hm = HashedMonomial(hexp)
    (haskey(ht.expmap, hm)) && (return hm)

    ht.expmap[hm] = exp
    hm
end


function add_product!(ht::MonomialHashtable, h1, h2)
    (haskey(ht.prods, (h1, h2))) && (return ht.prods[(h1, h2)])

    # make hash additive
    prode = ht.expmap[h1] .+ ht.expmap[h2]
    prodh = add_monomial!(ht, prode)

    ht.prods[(h1, h2)] = prodh

    prodh
end


# returns exponent vector for the given hash
function Base.getindex(ht::MonomialHashtable, h::HashedMonomial)
    getindex(ht.expmap, h)
end


#------------------------------------------------------------------------------
