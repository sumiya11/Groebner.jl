# Pairset and Basis

# Basis is a structure that stores a list of polynomials. Each polynomial is
# represented with a sorted vector of monomials and a vector of coefficients.
# Monomials and coefficients are stored in the basis separately. Each monomial
# is represented with an integer --- an index to a bucket in the hashtable (see
# f4/hashtable.jl).

# Pairset is a list of critical pairs (SPairs).

# S-pair{Degree}, or a pair of polynomials,
struct SPair{Degree}
    # First polynomial given by its index in the basis array
    poly1::Int32
    # Second polynomial -//-
    poly2::Int32
    # Index of lcm(lead(poly1), lead(poly2)) in the hashtable
    lcm::MonomIdx
    # Total degree of lcm
    deg::Degree
end

# Stores S-Pairs and some additional info.
mutable struct Pairset{Degree}
    pairs::Vector{SPair{Degree}}
    # A buffer of monomials represented with indices to a hashtable
    lcms::Vector{MonomIdx}
    # Number of filled pairs, initially zero
    load::Int
end

# Initializes and returns a pairset with max_vars_in_monom for `initial_size` pairs.
function initialize_pairset(::Type{Degree}; initial_size=2^6) where {Degree}
    pairs = Vector{SPair{Degree}}(undef, initial_size)
    lcms = Vector{MonomIdx}(undef, 0)
    Pairset(pairs, lcms, 0)
end

Base.isempty(ps::Pairset) = ps.load == 0

# Checks if it is possible to add `to_add` number of pairs to the pairset, and
# resizes the pairset if not
function resize_pairset_if_needed!(ps::Pairset, to_add::Int)
    newsize = length(ps.pairs)
    while ps.load + to_add > newsize
        newsize = max(2 * newsize, ps.load + to_add)
    end
    resize!(ps.pairs, newsize)
    nothing
end

function resize_lcms_if_needed!(ps::Pairset, nfilled::Int)
    if length(ps.lcms) < nfilled + 1
        resize!(ps.lcms, floor(Int, nfilled * 1.1) + 1)
    end
    nothing
end

# The type of the sugar degree
const SugarCube = Int

# Stores basis generators and some additional info
mutable struct Basis{C <: Coeff}
    # Vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with its index in hashtable
    monoms::Vector{Vector{MonomIdx}}
    # Polynomial coefficients
    coeffs::Vector{Vector{C}}

    # Max. number of polynomials that the basis can hold
    size::Int
    # Number of processed polynomials, initially zero
    nprocessed::Int
    # Total number of polys filled, initially zero
    nfilled::Int

    # If element of the basis at some index is redundant
    isredundant::Vector{Bool}
    # Positions of non-redundant elements in the basis
    nonredundant::Vector{Int}
    # Division masks of leading monomials of non-redundant basis elements
    divmasks::Vector{DivisionMask}
    # The number of non redundant elements in the basis
    nnonredundant::Int

    # Sugar degrees of basis polynomials
    sugar_cubes::Vector{SugarCube}
end

# Initialize basis for `ngens` elements with coefficient of type T
function initialize_basis(ring::PolyRing, ngens::Int, ::Type{T}) where {T <: Coeff}
    sz = ngens # * 2
    ndone = 0
    nfilled = 0
    nlead = 0

    monoms = Vector{Vector{MonomIdx}}(undef, sz)
    coeffs = Vector{Vector{T}}(undef, sz)
    isred = zeros(Bool, sz)
    nonred = Vector{Int}(undef, sz)
    lead = Vector{DivisionMask}(undef, sz)
    sugar_cubes = Vector{Int}(undef, sz)

    Basis(monoms, coeffs, sz, ndone, nfilled, isred, nonred, lead, nlead, sugar_cubes)
end

# initialize basis with the given hashed monomials and coefficients.
function initialize_basis(
    ring::PolyRing,
    hashedexps::Vector{Vector{MonomIdx}},
    coeffs::Vector{Vector{T}}
) where {T <: Coeff}
    sz = length(hashedexps)
    ndone = 0
    nfilled = sz
    nlead = 0

    isred = zeros(Bool, sz)
    nonred = Vector{Int}(undef, sz)
    lead = Vector{DivisionMask}(undef, sz)
    sugar_cubes = Vector{Int}(undef, sz)

    Basis(hashedexps, coeffs, sz, ndone, nfilled, isred, nonred, lead, nlead, sugar_cubes)
end

function repr_basis(basis::Basis{T}) where {T}
    s = """
    $(typeof(basis))
    Filled / Processed / Non-redundant : $(basis.nfilled) / $(basis.nprocessed) / $(basis.nnonredundant)
    Size allocated: $(basis.size)"""
    s
end

function copy_basis(
    basis::Basis{C},
    new_coeffs::Vector{Vector{T}};
    deepcopy=true
) where {C <: Coeff, T <: Coeff}
    if deepcopy
        basis = _deepcopy_basis(basis, new_coeffs)
    end
    monoms = basis.monoms
    coeffs = new_coeffs
    isred = basis.isredundant
    nonred = basis.nonredundant
    divmasks = basis.divmasks
    sugar_cubes = basis.sugar_cubes
    Basis(
        monoms,
        coeffs,
        basis.size,
        basis.nprocessed,
        basis.nfilled,
        isred,
        nonred,
        divmasks,
        basis.nnonredundant,
        sugar_cubes
    )
end

function _deepcopy_basis(basis::Basis{T}, new_coeffs::Vector{Vector{C}}) where {T, C}
    monoms = Vector{Vector{MonomIdx}}(undef, basis.size)
    @inbounds for i in 1:(basis.nfilled)
        monoms[i] = Vector{MonomIdx}(undef, length(basis.monoms[i]))
        for j in 1:length(basis.monoms[i])
            monoms[i][j] = basis.monoms[i][j]
        end
    end
    isred = copy(basis.isredundant)
    nonred = copy(basis.nonredundant)
    lead = copy(basis.divmasks)
    sugar_cubes = copy(basis.sugar_cubes)
    Basis(
        monoms,
        new_coeffs,
        basis.size,
        basis.nprocessed,
        basis.nfilled,
        isred,
        nonred,
        lead,
        basis.nnonredundant,
        sugar_cubes
    )
end

function deepcopy_basis(basis::Basis{T}) where {T}
    new_coeffs = Vector{Vector{T}}(undef, basis.size)
    @inbounds for i in 1:(basis.nfilled)
        new_coeffs[i] = Vector{T}(undef, length(basis.coeffs[i]))
        for j in 1:length(basis.coeffs[i])
            new_coeffs[i][j] = copy(basis.coeffs[i][j])
        end
    end
    _deepcopy_basis(basis, new_coeffs)
end

# 
function resize_basis_if_needed!(basis::Basis{T}, to_add::Int) where {T}
    while basis.nprocessed + to_add >= basis.size
        basis.size = max(basis.size * 2, basis.nprocessed + to_add)
        resize!(basis.monoms, basis.size)
        resize!(basis.coeffs, basis.size)
        resize!(basis.isredundant, basis.size)
        @inbounds basis.isredundant[(basis.nprocessed + 1):end] .= false
        resize!(basis.nonredundant, basis.size)
        resize!(basis.divmasks, basis.size)
        resize!(basis.sugar_cubes, basis.size)
    end
    @invariant basis.size >= basis.nprocessed + to_add
    nothing
end

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(ring, basis::Basis{<:CoeffFF})
    @log level = -4 "Normalizing polynomials in the basis"
    cfs = basis.coeffs
    @inbounds for i in 1:(basis.nfilled)
        !isassigned(cfs, i) && continue   # TODO: this is kind of bad
        ch = ring.ch
        mul = invmod(cfs[i][1], ch) % ch
        for j in 2:length(cfs[i])
            cfs[i][j] = (cfs[i][j] * mul) % ch
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

# Normalize each element of the input basis
# by dividing it by leading coefficient
function normalize_basis!(ring, basis::Basis{<:CoeffQQ})
    @log level = -4 "Normalizing polynomials in the basis"
    cfs = basis.coeffs
    @inbounds for i in 1:(basis.nfilled)
        !isassigned(cfs, i) && continue
        mul = inv(cfs[i][1])
        for j in 2:length(cfs[i])
            cfs[i][j] *= mul
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

# Generate new S-pairs from pairs of polynomials
#   (basis[idx], basis[i])
# for every i < idx
#
# NOTE: discarding redundant critical pairs is unaffected by the current
# selection strategy
function update_pairset!(
    pairset::Pairset,
    basis::Basis{C},
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {C <: Coeff, M <: Monom}
    pr = entrytype(M)
    pl, bl = pairset.load, idx
    ps = pairset.pairs
    lcms = pairset.lcms

    new_lead = basis.monoms[idx][1]

    # generate a pair for each pair
    @inbounds for i in 1:(bl - 1)
        newidx = pl + i
        if !basis.isredundant[i] &&
           !is_gcd_const(ht.monoms[basis.monoms[i][1]], ht.monoms[new_lead])
            lcms[i] = get_lcm(basis.monoms[i][1], new_lead, ht, update_ht)
            deg = update_ht.hashdata[lcms[i]].deg
            ps[newidx] = SPair(Int32(i), Int32(idx), lcms[i], pr(deg))
        else
            # lcm == 0 will mark redundancy of an S-pair
            lcms[i] = MonomIdx(0)
            # ps[newidx] = SPair(i, idx, MonomIdx(0), pr(deg))
            ps[newidx] = SPair{pr}(Int32(i), Int32(idx), MonomIdx(0), typemax(pr))
        end
    end

    # traverse existing pairs
    @inbounds for i in 1:pl
        if iszero(ps[i].lcm)
            continue
        end

        j = ps[i].poly1
        l = ps[i].poly2
        m = max(ps[pl + l].deg, ps[pl + j].deg)

        # if an existing pair is divisible by the lead of new poly
        # and has a greater degree than newly generated one then
        if ps[i].deg > m && is_monom_divisible(ps[i].lcm, new_lead, ht)
            # mark an existing pair redundant
            ps[i] = SPair{pr}(ps[i].poly1, ps[i].poly2, MonomIdx(0), ps[i].deg)
        end
    end

    # traverse new pairs to move not-redundant ones first 
    j = 1
    @inbounds for i in 1:(bl - 1)
        if !basis.isredundant[i]
            ps[pl + j] = ps[pl + i]
            j += 1
        end
    end

    sort_pairset_by_degree!(pairset, pl + 1, j - 2)

    @inbounds for i in 1:(j - 1)
        lcms[i] = ps[pl + i].lcm
    end
    @inbounds lcms[j] = 0
    pc = j
    pc -= 1

    # mark redundancy of some pairs from lcms array
    @inbounds for j in 1:pc
        # if is not redundant already
        if !iszero(lcms[j])
            check_monomial_division_in_update(lcms, j + 1, pc, lcms[j], update_ht)
        end
    end

    # remove useless pairs from pairset
    # by moving them to the end
    j = 1
    @inbounds for i in 1:(pairset.load)
        iszero(ps[i].lcm) && continue
        ps[j] = ps[i]
        j += 1
    end

    # ensure that basis hashtable can store new lcms
    resize_hashtable_if_needed!(ht, pc)

    # add new lcms to the basis hashtable,
    # including index j and not including index pc
    insert_lcms_in_basis_hash_table!(pairset, pl, ht, update_ht, basis, lcms, j, pc + 1)

    # mark redundant polynomials in basis
    nonred = basis.nonredundant
    lml = basis.nnonredundant
    @inbounds for i in 1:lml
        if !basis.isredundant[nonred[i]]
            if is_monom_divisible(basis.monoms[nonred[i]][1], new_lead, ht)
                basis.isredundant[nonred[i]] = true
            end
        end
    end

    nothing
end

# Updates information about redundant generators in the basis
function update_basis!(basis::Basis, ht::MonomialHashtable{M}) where {M <: Monom}
    k = 1
    lead = basis.divmasks
    nonred = basis.nonredundant
    @inbounds for i in 1:(basis.nnonredundant)
        if !basis.isredundant[nonred[i]]
            basis.divmasks[k] = lead[i]
            basis.nonredundant[k] = nonred[i]
            k += 1
        end
    end
    basis.nnonredundant = k - 1

    @inbounds for i in (basis.nprocessed + 1):(basis.nfilled)
        if !basis.isredundant[i]
            lead[k] = ht.hashdata[basis.monoms[i][1]].divmask
            nonred[k] = i
            k += 1
        end
    end

    basis.nnonredundant = k - 1
    basis.nprocessed = basis.nfilled
end

# Checks if element of basis at position idx is redundant
function is_redundant!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {M}
    pt = entrytype(M)
    resize_hashtable_if_needed!(update_ht, 0)

    # lead of new polynomial
    lead_new = basis.monoms[idx][1]
    ps = pairset.pairs

    @inbounds for i in (idx + 1):(basis.nfilled)
        i == idx && continue
        basis.isredundant[i] && continue

        # lead of new polynomial at index i > idx
        lead_i = basis.monoms[i][1]

        if is_monom_divisible(lead_new, lead_i, ht)
            # add new S-pair corresponding to Spoly(i, idx)
            lcm_new = get_lcm(lead_i, lead_new, ht, ht)
            psidx = pairset.load + 1
            ps[psidx] =
                SPair{pt}(Int32(i), Int32(idx), lcm_new, pt(ht.hashdata[lcm_new].deg))

            # mark redundant
            basis.isredundant[idx] = true
            pairset.load += 1

            return true
        end
    end

    return false
end

# Updates basis and pairset.
# 
# New elements added to the basis from the f4 matrix 
# (the ones with indices from basis.nprocessed+1 to basis.nfilled)
# are checked for being redundant.
#
# Then, pairset is updated with the new S-pairs formed
# from the added elements
#
# Always checks new elements redundancy.
@timed_block function update!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M}
) where {M <: Monom}
    # total number of elements in the basis (old + new)
    npivs = basis.nfilled
    # number of potential critical pairs to add
    npairs = basis.nprocessed * npivs + div((npivs + 1) * npivs, 2)

    # make sure pairset and update hashtable have enough
    # space to store new pairs
    # note: we create too big array, can be fixed
    resize_pairset_if_needed!(pairset, npairs)
    pairset_size = length(pairset.pairs)

    # update pairset,
    # for each new element in basis
    @inbounds for i in (basis.nprocessed + 1):(basis.nfilled)
        # check redundancy of new polynomial
        is_redundant!(pairset, basis, ht, update_ht, i) && continue
        resize_lcms_if_needed!(pairset, basis.nfilled)
        # if not redundant, then add new S-pairs to pairset
        update_pairset!(pairset, basis, ht, update_ht, i)
    end

    # update basis
    update_basis!(basis, ht)
    pairset_size
end

# given input exponent and coefficient vectors hashes exponents into `ht`
# and then constructs hashed polynomials for `basis`
function fill_data!(
    basis::Basis,
    ht::MonomialHashtable{M},
    exponents::Vector{Vector{M}},
    coeffs::Vector{Vector{T}}
) where {M, T}
    ngens = length(exponents)
    @inbounds for i in 1:ngens
        resize_hashtable_if_needed!(ht, length(exponents[i]))

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

    basis.nfilled = ngens
end

function sweep_redundant!(basis::Basis, hashtable)
    # here -- assert that basis is in fact a Groebner basis.
    # NOTE: maybe sort generators for more effective sweeping?
    @inbounds for i in 1:(basis.nprocessed)
        for j in (i + 1):(basis.nprocessed)
            basis.isredundant[i] && continue
            basis.isredundant[j] && continue
            lead_i = basis.monoms[i][1]
            lead_j = basis.monoms[j][1]
            if is_monom_divisible(lead_i, lead_j, hashtable)
                basis.isredundant[i] = true
            elseif is_monom_divisible(lead_j, lead_i, hashtable)
                basis.isredundant[j] = true
            end
        end
    end
    nothing
end

# Remove redundant elements from the basis
# by moving all non-redundant up front
function mark_redundant!(basis::Basis)
    j = 1
    @inbounds for i in 1:(basis.nnonredundant)
        if !basis.isredundant[basis.nonredundant[i]]
            basis.divmasks[j] = basis.divmasks[i]
            basis.nonredundant[j] = basis.nonredundant[i]
            j += 1
        end
    end
    basis.nnonredundant = j - 1
    @assert basis.nprocessed == basis.nfilled
    basis
end

# f4 must produce a `Basis` object which satisfies several conditions.
# This functions standardizes the given basis so that conditions hold.
# (see f4/f4.jl)
function standardize_basis!(ring, basis::Basis, ht::MonomialHashtable, ord)
    @inbounds for i in 1:(basis.nnonredundant)
        idx = basis.nonredundant[i]
        basis.nonredundant[i] = i
        basis.isredundant[i] = false
        basis.coeffs[i] = basis.coeffs[idx]
        basis.monoms[i] = basis.monoms[idx]
    end
    basis.size = basis.nprocessed = basis.nfilled = basis.nnonredundant
    resize!(basis.coeffs, basis.nprocessed)
    resize!(basis.monoms, basis.nprocessed)
    resize!(basis.divmasks, basis.nprocessed)
    resize!(basis.nonredundant, basis.nprocessed)
    resize!(basis.isredundant, basis.nprocessed)
    resize!(basis.sugar_cubes, basis.nprocessed)
    sort_polys_by_lead_increasing!(basis, ht, ord=ord)
    normalize_basis!(ring, basis)
end

# Returns the exponent vectors of polynomials in the basis
function hash_to_exponents(basis::Basis, ht::MonomialHashtable{M}) where {M <: Monom}
    exps = Vector{Vector{M}}(undef, basis.nnonredundant)
    @inbounds for i in 1:(basis.nnonredundant)
        idx = basis.nonredundant[i]
        poly = basis.monoms[idx]
        exps[i] = Vector{M}(undef, length(poly))
        for j in 1:length(poly)
            exps[i][j] = ht.monoms[poly[j]]
        end
    end
    exps
end

# Returns the exponent vectors and the coefficients of polynomials in the basis
function export_basis_data(
    basis::Basis{C},
    ht::MonomialHashtable{M}
) where {M <: Monom, C <: Coeff}
    exps = hash_to_exponents(basis, ht)
    coeffs = Vector{Vector{C}}(undef, basis.nnonredundant)
    @inbounds for i in 1:(basis.nnonredundant)
        idx = basis.nonredundant[i]
        coeffs[i] = basis.coeffs[idx]
    end
    exps, coeffs
end

# For a given list of S-pairs and a list of indices `plcm`
# adds indices from plcm[ifirst:ilast]
# to the hashtable ht
function insert_lcms_in_basis_hash_table!(
    pairset::Pairset,
    off::Int,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    basis::Basis,
    plcm::Vector{MonomIdx},
    ifirst::Int,
    ilast::Int
) where {M}
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

        if is_gcd_const(monoms[ps[off + l].poly1][1], monoms[ps[off + 1].poly2][1], ht)
            l += 1
            continue
        end

        ps[m] = ps[off + l]

        h = update_ht.hashdata[plcm[l]].hash
        ht.monoms[ht.load + 1] = copy_monom(update_ht.monoms[plcm[l]])
        n = ht.monoms[ht.load + 1]

        k = h
        i = MonomHash(1)
        @inbounds while i <= ht.size
            k = next_lookup_index(h, i, mod)
            hm = ht.hashtable[k]

            # if free
            iszero(hm) && break

            if ishashcollision(ht, hm, n, h)
                i += MonomHash(1)
                continue
            end

            ps[m] = SPair{typeof(ps[m].deg)}(ps[m].poly1, ps[m].poly2, hm, ps[m].deg)
            m += 1
            l += 1
            @goto Letsgo
        end

        ht.hashtable[k] = pos = ht.load + 1

        uhd = update_ht.hashdata
        ll = plcm[l]
        ht.hashdata[ht.load + 1] = Hashvalue(0, h, uhd[ll].divmask, uhd[ll].deg)

        ht.load += 1
        ps[m] = SPair{typeof(ps[m].deg)}(ps[m].poly1, ps[m].poly2, MonomIdx(pos), ps[m].deg)
        m += 1
        l += 1
    end

    pairset.load = m - 1
end
