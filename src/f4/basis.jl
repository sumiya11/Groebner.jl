# Pairset and Basis

# Basis is a structure that stores a list of polynomials in f4.
# Each polynomial is represented as a sorted vector of monomials 
# and a sorted vector of coefficients.
# Monomials and coefficients are stored in the basis separately. 
# Each monomial is represented with an integer -- 
# an index to a bucket in the hashtable (see f4/hashtable.jl)

# Pairset is a list of SPairs.
# An SPair is a critical pair in f4.

# S-pair{Power}, a pair of polynomials,
# Power is the type of exponent vector entry
struct SPair{Power<:Integer}
    # first polynomial as an index from the basis array
    poly1::Int
    # second polynomial -//-
    poly2::Int
    # index of lcm(lt(poly1), lt(poly2)) in hashtable
    lcm::MonomIdx
    # total degree of lcm
    deg::Power
end

# Stores S-Pairs
mutable struct Pairset{Power}
    pairs::Vector{SPair{Power}}
    # number of filled pairs,
    # Initially zero
    load::Int
end

# Initializes and returns pairset with capacity for `initial_size` pairs.
# fine tuned parameter 2^6
function initialize_pairset(::Type{Power}; initial_size=2^6) where {Power} 
    pairs = Vector{SPair{Power}}(undef, initial_size)
    Pairset(pairs, 0)
end

Base.isempty(ps::Pairset) = ps.load == 0

# Checks if it is possible to add `added` number of pairs to the pairset,
# and extends the pairset if not
function check_enlarge_pairset!(ps::Pairset, added::Int)
    sz = length(ps.pairs)
    if ps.load + added >= sz
        newsz = max(2 * sz, ps.load + added)
        resize!(ps.pairs, newsz)
    end
    nothing
end

#------------------------------------------------------------------------------

# Stores basis generators and some additional info
mutable struct Basis{C<:Coeff}
    # vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with its index in hashtable
    monoms::Vector{Vector{MonomIdx}}
    # polynomial coefficients
    coeffs::Vector{Vector{C}}

    #= Keeping track of sizes   =#
    #=  ndone <= ntotal <= size =#
    # total allocated size
    size::Int
    # number of processed polynomials in `monoms`
    # Iitially zero
    ndone::Int
    # total number of polys filled in `monoms`
    # (these will be handled during next update! call)
    # Iitially zero
    ntotal::Int

    # always true
    # length(monoms) == length(coeffs) == length(isred) == size

    #= Keeping track of redundancy =#
    #= invariant =#
    #= length(lead) == length(nonred) == count(isred) == nlead =#
    # if element of the basis
    # is redundant
    isred::Vector{Bool}
    # positions of non-redundant elements in the basis
    nonred::Vector{Int}
    # division masks of leading monomials of
    # non redundant basis elements
    lead::Vector{DivisionMask}
    # number of filled elements in lead
    nlead::Int
end

# initialize basis for `ngens` elements with coefficient of type T
function initialize_basis(ring::PolyRing, ngens::Int, ::Type{T}) where {T<:Coeff}
    sz = ngens * 2
    ndone = 0
    ntotal = 0
    nlead = 0

    monoms = Vector{Vector{MonomIdx}}(undef, sz)
    coeffs = Vector{Vector{T}}(undef, sz)
    isred = zeros(Bool, sz)
    nonred = Vector{Int}(undef, sz)
    lead = Vector{DivisionMask}(undef, sz)

    Basis(monoms, coeffs, sz, ndone, ntotal, isred, nonred, lead, nlead)
end

# initialize basis with the given hashed monomials and coefficients
function initialize_basis(ring::PolyRing, hashedexps::Vector{Vector{MonomIdx}}, coeffs::Vector{Vector{T}}) where {T<:Coeff}
    sz = length(hashedexps)
    ndone = 0
    ntotal = 0
    nlead = 0

    isred = zeros(Bool, sz)
    nonred = Vector{Int}(undef, sz)
    lead = Vector{DivisionMask}(undef, sz)

    Basis(hashedexps, coeffs, sz, ndone, ntotal, isred, nonred, lead, nlead)
end

function deepcopy_basis(basis::Basis{T}) where {T}
    # why not extend Base.deepcopy ?...
    monoms = Vector{Vector{MonomIdx}}(undef, basis.size)
    coeffs = Vector{Vector{T}}(undef, basis.size)
    @inbounds for i in 1:basis.ntotal
        monoms[i] = Vector{MonomIdx}(undef, length(basis.monoms[i]))
        coeffs[i] = Vector{T}(undef, length(basis.coeffs[i]))
        @inbounds for j in 1:length(basis.monoms[i])
            monoms[i][j] = basis.monoms[i][j]
            coeffs[i][j] = basis.coeffs[i][j]
        end
    end
    isred = copy(basis.isred)
    nonred = copy(basis.nonred)
    lead = copy(basis.lead)
    Basis(monoms, coeffs, basis.size, basis.ndone,
        basis.ntotal, isred, nonred, lead,
        basis.nlead)
end

function check_enlarge_basis!(basis::Basis{T}, added::Int) where {T}
    while basis.ndone + added >= basis.size
        basis.size = max(basis.size * 2, basis.ndone + added)
        resize!(basis.monoms, basis.size)
        resize!(basis.coeffs, basis.size)
        resize!(basis.isred, basis.size)
        @inbounds basis.isred[basis.ndone+1:end] .= false
        resize!(basis.nonred, basis.size)
        resize!(basis.lead, basis.size)
    end
    nothing
end

#------------------------------------------------------------------------------

function cleanup_basis!(ring::PolyRing, basis::Basis, prime)
    ring.ch = prime
    normalize_basis!(ring, basis)
end

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(ring, basis::Basis{<:CoeffFF})
    cfs = basis.coeffs
    @inbounds for i in 1:basis.ntotal
        !isassigned(cfs, i) && continue
        ch = ring.ch
        mul = invmod(cfs[i][1], ch) % ch
        @inbounds for j in 2:length(cfs[i])
            cfs[i][j] = (cfs[i][j] * mul) % ch
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(ring, basis::Basis{<:CoeffQQ})
    cfs = basis.coeffs
    @inbounds for i in 1:basis.ntotal
        !isassigned(cfs, i) && continue
        mul = inv(cfs[i][1])
        @inbounds for j in 2:length(cfs[i])
            cfs[i][j] *= mul
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

#------------------------------------------------------------------------------

# Generate new S-pairs from pairs of polynomials
#   (basis[idx], basis[i])
# for every i < idx
function update_pairset!(
        pairset::Pairset,
        basis::Basis{C},
        ht::MonomialHashtable{M},
        update_ht::MonomialHashtable{M},
        idx::Int,
        plcm::Vector{MonomIdx}) where {C<:Coeff, M<:Monom}

    pr = powertype(M)
    pl, bl = pairset.load, idx
    ps = pairset.pairs

    new_lead = basis.monoms[idx][1]

    # generate a pair for each pair
    @inbounds for i in 1:bl-1
        plcm[i] = get_lcm(basis.monoms[i][1], new_lead, ht, update_ht)
        deg = update_ht.hashdata[plcm[i]].deg
        newidx = pl + i
        if !basis.isred[i]
            ps[newidx] = SPair(i, idx, plcm[i], pr(deg))
        else
            # lcm == 0 will mark redundancy of an S-pair
            ps[newidx] = SPair(i, idx, MonomIdx(0), pr(deg))
        end
    end

    # traverse existing pairs
    @inbounds for i in 1:pl
        j = ps[i].poly1
        l = ps[i].poly2
        m = max(ps[pl+l].deg, ps[pl+j].deg)

        # if an existing pair is divisible by the lead of new poly
        # and has a greater degree than newly generated one then
        if is_monom_divisible(ps[i].lcm, new_lead, ht) && ps[i].deg > m
            # mark an existing pair redundant
            ps[i] = SPair(ps[i].poly1, ps[i].poly2, MonomIdx(0), ps[i].deg)
        end
    end

    # traverse new pairs to move not-redundant ones first 
    j = 1
    @inbounds for i in 1:bl-1
        if !basis.isred[i]
            ps[pl+j] = ps[pl+i]
            j += 1
        end
    end

    sort_pairset_by_degree!(pairset, pl + 1, j - 2)

    @inbounds for i in 1:j-1
        plcm[i] = ps[pl+i].lcm
    end
    @inbounds plcm[j] = 0
    pc = j
    pc -= 1

    # mark redundancy of some pairs from plcm array
    @inbounds for j in 1:pc
        # if is not redundant already
        if !iszero(plcm[j])
            check_monomial_division_in_update(plcm, j + 1, pc, plcm[j], update_ht)
        end
    end

    # remove useless pairs from pairset
    # by moving them to the end
    j = 1
    @inbounds for i in 1:pairset.load
        iszero(ps[i].lcm) && continue
        ps[j] = ps[i]
        j += 1
    end

    # ensure that basis hashtable can store new lcms
    # YES
    check_enlarge_hashtable!(ht, pc)
    # if ht.size - ht.load <= pc
    #     enlarge_hash_table!(ht)
    # end

    # add new lcms to the basis hashtable,
    # including index j and not including index pc
    insert_plcms_in_basis_hash_table!(pairset, pl, ht, update_ht, basis, plcm, j, pc + 1)

    # mark redundant polynomials in basis
    nonred = basis.nonred
    lml = basis.nlead
    @inbounds for i in 1:lml
        if !basis.isred[nonred[i]]
            if is_monom_divisible(basis.monoms[nonred[i]][1], new_lead, ht)
                basis.isred[nonred[i]] = true
            end
        end
    end

    nothing
end

# Updates information about redundant generators in the basis
function update_basis!(
        basis::Basis,
        ht::MonomialHashtable{M}) where {M<:Monom}

    k = 1
    lead = basis.lead
    nonred = basis.nonred
    @inbounds for i in 1:basis.nlead
        if !basis.isred[nonred[i]]
            basis.lead[k] = lead[i]
            basis.nonred[k] = nonred[i]
            k += 1
        end
    end
    basis.nlead = k - 1

    @inbounds for i in basis.ndone+1:basis.ntotal
        if !basis.isred[i]
            lead[k] = ht.hashdata[basis.monoms[i][1]].divmask
            nonred[k] = i
            k += 1
        end
    end

    basis.nlead = k - 1
    basis.ndone = basis.ntotal
end

# Checks if element of basis at position idx is redundant
function is_redundant!(
        pairset::Pairset, basis::Basis, ht::MonomialHashtable{M}, 
        update_ht::MonomialHashtable{M}, idx::Int) where {M}

    pt = powertype(M)
    check_enlarge_hashtable!(update_ht, 0)

    # lead of new polynomial
    lead_new = basis.monoms[idx][1]
    ps = pairset.pairs

    @inbounds for i in idx+1:basis.ntotal
        i == idx && continue
        basis.isred[i] && continue

        # lead of new polynomial at index i > idx
        lead_i = basis.monoms[i][1]

        if is_monom_divisible(lead_new, lead_i, ht)
            # add new S-pair corresponding to Spoly(i, idx)
            lcm_new = get_lcm(lead_i, lead_new, ht, ht)
            psidx = pairset.load + 1
            ps[psidx] = SPair(i, idx, lcm_new, pt(ht.hashdata[lcm_new].deg))

            # mark redundant
            basis.isred[idx] = true
            pairset.load += 1

            return true
        end
    end

    return false
end

function check_enlarge_plcm!(plcm::Vector{T}, ntotal) where {T}
    if length(plcm) < ntotal + 1
        resize!(plcm, floor(Int, ntotal * 1.1) + 1)
    end
    nothing
end

# Updates basis and pairset.
# 
# New elements added to the basis from the f4 matrix 
# (the ones with indices from basis.ndone+1 to basis.ntotal)
# are checked for being redundant.
#
# Then, pairset is updated with the new S-pairs formed
# from the added elements
#
# Always checks new elements redundancy.
function update!(
        pairset::Pairset,
        basis::Basis,
        ht::MonomialHashtable{M},
        update_ht::MonomialHashtable{M},
        plcm::Vector{MonomIdx}) where {M<:Monom}

    # total number of elements in the basis (old + new)
    npivs = basis.ntotal
    # number of potential critical pairs to add
    npairs = basis.ndone * npivs + div((npivs + 1)*npivs, 2)

    # make sure pairset and update hashtable have enough
    # space to store new pairs
    # note: we create too big array, can be fixed
    check_enlarge_pairset!(pairset, npairs)
    pairset_size = length(pairset.pairs)

    # update pairset,
    # for each new element in basis
    @inbounds for i in basis.ndone+1:basis.ntotal
        # check redundancy of new polynomial
        is_redundant!(pairset, basis, ht, update_ht, i) && continue
        # YES
        check_enlarge_plcm!(plcm, basis.ntotal)
        # if not redundant, then add new S-pairs to pairset
        update_pairset!(pairset, basis, ht, update_ht, i, plcm)
    end

    # update basis
    update_basis!(basis, ht)

    pairset_size
end

#------------------------------------------------------------------------------

# given input exponent and coefficient vectors hashes exponents into `ht`
# and then constructs hashed polynomials for `basis`
function fill_data!(basis::Basis, ht::MonomialHashtable{M}, 
        exponents::Vector{Vector{M}}, coeffs::Vector{Vector{T}}) where {M,T}
    
    ngens = length(exponents)
    @inbounds for i in 1:ngens
        check_enlarge_hashtable!(ht, length(exponents[i]))

        nterms = length(coeffs[i])
        basis.coeffs[i] = coeffs[i]
        basis.monoms[i] = Vector{MonomIdx}(undef, nterms)
        poly = basis.monoms[i]
        @inbounds for j in 1:nterms
            poly[j] = insert_in_hash_table!(ht, exponents[i][j])
        end

        # sort terms (not needed)
        # beautifuly coefficients (not needed)
    end

    basis.ntotal = ngens
end

# Remove redundant elements from the basis
# by moving all non-redundant up front
function filter_redundant!(basis::Basis)
    j = 1
    @inbounds for i in 1:basis.nlead
        if !basis.isred[basis.nonred[i]]
            basis.lead[j] = basis.lead[i]
            basis.nonred[j] = basis.nonred[i]
            j += 1
        end
    end
    basis.nlead = j - 1
    @assert basis.ndone == basis.ntotal
    basis
end

# f4 must produce a `Basis` object which satisfies several conditions.
# This functions standardizes the given basis so that conditions hold.
# (see f4/f4.jl)
function standardize_basis!(ring, basis::Basis, ht::MonomialHashtable, ord)
    @inbounds for i in 1:basis.nlead
        idx = basis.nonred[i]
        basis.nonred[i] = i
        basis.isred[i] = false
        basis.coeffs[i] = basis.coeffs[idx]
        basis.monoms[i] = basis.monoms[idx]
    end
    basis.size = basis.ndone = basis.ntotal = basis.nlead
    resize!(basis.coeffs, basis.ndone)
    resize!(basis.monoms, basis.ndone)
    resize!(basis.lead, basis.ndone)
    resize!(basis.nonred, basis.ndone)
    resize!(basis.isred, basis.ndone)
    sort_gens_by_lead_increasing!(basis, ht, ord=ord)
    normalize_basis!(ring, basis)
end

# Returns the exponent vectors of polynomials in the basis
function hash_to_exponents(basis::Basis, ht::MonomialHashtable{M}) where {M<:Monom}
    exps = Vector{Vector{M}}(undef, basis.nlead)
    @inbounds for i in 1:basis.nlead
        idx = basis.nonred[i]
        poly = basis.monoms[idx]
        exps[i] = Vector{M}(undef, length(poly))
        @inbounds for j in 1:length(poly)
            exps[i][j] = ht.exponents[poly[j]]
        end
    end
    exps
end

# Returns the exponent vectors and the coefficients of polynomials in the basis
function export_basis_data(basis::Basis{C}, ht::MonomialHashtable{M}) where {M<:Monom,C<:Coeff}
    exps = hash_to_exponents(basis, ht)
    coeffs = Vector{Vector{C}}(undef, basis.nlead)
    @inbounds for i in 1:basis.nlead
        idx = basis.nonred[i]
        coeffs[i] = basis.coeffs[idx]
    end
    exps, coeffs
end

#------------------------------------------------------------------------------

# For a given list of S-pairs and a list of indices `plcm`
# adds indices from plcm[ifirst:ilast]
# to the hashtable ht
function insert_plcms_in_basis_hash_table!(
        pairset::Pairset,
        off::Int,
        ht::MonomialHashtable{M},
        update_ht::MonomialHashtable{M},
        basis::Basis,
        plcm::Vector{MonomIdx},
        ifirst::Int, ilast::Int) where {M}
    # including ifirst and not including ilast

    monoms = basis.monoms
    ps = pairset.pairs

    mod = MonomHash(ht.size - 1)
    @assert ispow2(mod + 1)

    m = ifirst
    l = 1
    @label Letsgo
    @inbounds while l < ilast
        if iszero(plcm[l])
            l += 1
            continue
        end

        if is_gcd_const(monoms[ps[off+l].poly1][1], monoms[ps[off+1].poly2][1], ht)
            l += 1
            continue
        end

        ps[m] = ps[off+l]

        h = update_ht.hashdata[plcm[l]].hash
        ht.exponents[ht.load+1] = copy(update_ht.exponents[plcm[l]])
        n = ht.exponents[ht.load+1]

        k = h
        i = MonomHash(1)
        @inbounds while i <= ht.size
            k = nexthashindex(h, i, mod)
            hm = ht.hashtable[k]
            
            # if free
            iszero(hm) && break

            if ishashcollision(ht, hm, n, h)
                i += MonomHash(1)
                continue
            end

            ps[m] = SPair(ps[m].poly1, ps[m].poly2, hm, ps[m].deg)
            m += 1
            l += 1
            @goto Letsgo
        end

        ht.hashtable[k] = pos = ht.load + 1

        uhd = update_ht.hashdata
        ll = plcm[l]
        ht.hashdata[ht.load+1] = Hashvalue(0, h, uhd[ll].divmask, uhd[ll].deg)

        ht.load += 1
        ps[m] = SPair(ps[m].poly1, ps[m].poly2, MonomIdx(pos), ps[m].deg)
        m += 1
        l += 1
    end

    pairset.load = m - 1
end
