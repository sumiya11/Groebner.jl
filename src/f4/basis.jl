# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# Pairset and Basis

# Pairset is a list of critical pairs (SPairs).

# Basis is a structure that stores a list of polynomials. Each polynomial is
# represented with a sorted vector of monomials and a vector of coefficients.
# Monomials and coefficients are stored in the basis separately. Each monomial
# is represented with an integer -- a unique identifier that indexes a bucket in
# the hashtable (see f4/hashtable.jl).

###
# Pairset

# S-pair{Degree}, or, a pair of polynomials,
struct CriticalPair{Degree}
    # First polynomial given by its index in the basis array
    poly1::Int32
    # Second polynomial -//-
    poly2::Int32
    # Index of lcm(lead(poly1), lead(poly2)) in the hashtable
    lcm::MonomId
    # Total degree of lcm
    deg::Degree
end

# Stores S-Pairs and some additional info.
mutable struct Pairset{Degree}
    pairs::Vector{CriticalPair{Degree}}
    # A buffer of monomials represented with indices to a hashtable
    lcms::Vector{MonomId}
    # Number of filled pairs, initially zero
    load::Int
end

# Initializes and returns a pairset with monom_max_vars for `initial_size` pairs.
function pairset_initialize(::Type{Degree}; initial_size=2^6) where {Degree}
    pairs = Vector{CriticalPair{Degree}}(undef, initial_size)
    lcms = Vector{MonomId}(undef, 0)
    Pairset(pairs, lcms, 0)
end

Base.isempty(ps::Pairset) = ps.load == 0

# Checks if it is possible to add `to_add` number of pairs to the pairset, and
# resizes the pairset if not
function pairset_resize_if_needed!(ps::Pairset, to_add::Int)
    newsize = length(ps.pairs)
    while ps.load + to_add > newsize
        newsize = max(2 * newsize, ps.load + to_add)
    end
    resize!(ps.pairs, newsize)
    nothing
end

function pairset_resize_lcms_if_needed!(ps::Pairset, nfilled::Int)
    if length(ps.lcms) < nfilled + 1
        resize!(ps.lcms, floor(Int, nfilled * 1.1) + 1)
    end
    nothing
end

###
# Basis

# The type of the sugar degree
const SugarCube = Int

# Stores basis generators and some additional info
mutable struct Basis{C <: Coeff}
    # Vector of polynomials, each polynomial is a vector of monomials,
    # each monomial is represented with its index in hashtable
    monoms::Vector{Vector{MonomId}}
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
function basis_initialize(ring::PolyRing, ngens::Int, ::Type{T}) where {T <: Coeff}
    sz = ngens # * 2
    ndone = 0
    nfilled = 0
    nlead = 0

    monoms = Vector{Vector{MonomId}}(undef, sz)
    coeffs = Vector{Vector{T}}(undef, sz)
    isred = zeros(Bool, sz)
    nonred = Vector{Int}(undef, sz)
    lead = Vector{DivisionMask}(undef, sz)
    sugar_cubes = Vector{Int}(undef, sz)

    Basis(monoms, coeffs, sz, ndone, nfilled, isred, nonred, lead, nlead, sugar_cubes)
end

# initialize basis with the given (already hashed) monomials and coefficients.
function basis_initialize(
    ring::PolyRing,
    hashedexps::Vector{Vector{MonomId}},
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

# Same as f4_initialize_structs, but uses an existing hashtable
function basis_initialize_using_existing_hashtable(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable;
) where {M, C <: Coeff}
    basis = basis_initialize(ring, length(monoms), C)
    fill_data!(basis, present_ht, monoms, coeffs)
    basis
end

###
# Basis utils

function basis_shallow_copy_with_new_coeffs(
    basis::Basis{C},
    new_coeffs::Vector{Vector{T}}
) where {C <: Coeff, T <: Coeff}
    Basis(
        basis.monoms,
        new_coeffs,
        basis.size,
        basis.nprocessed,
        basis.nfilled,
        basis.isredundant,
        basis.nonredundant,
        basis.divmasks,
        basis.nnonredundant,
        basis.sugar_cubes
    )
end

function basis_deep_copy_with_new_coeffs(
    basis::Basis{C},
    new_coeffs::Vector{Vector{T}}
) where {C <: Coeff, T <: Coeff}
    # Assume that MonomId is trivially copiable    
    @invariant isbitstype(MonomId)

    monoms = Vector{Vector{MonomId}}(undef, length(basis.monoms))
    @inbounds for i in 1:length(basis.monoms)
        !isassigned(basis.monoms, i) && continue
        monoms[i] = Vector{MonomId}(undef, length(basis.monoms[i]))
        for j in 1:length(basis.monoms[i])
            monoms[i][j] = basis.monoms[i][j]
        end
    end

    Basis(
        monoms,
        new_coeffs,
        basis.size,
        basis.nprocessed,
        basis.nfilled,
        copy(basis.isredundant),
        copy(basis.nonredundant),
        copy(basis.divmasks),
        basis.nnonredundant,
        copy(basis.sugar_cubes)
    )
end

function basis_deepcopy(basis::Basis{C}) where {C <: Coeff}
    coeffs = Vector{Vector{C}}(undef, length(basis.coeffs))

    if isbitstype(C) # For Z/pZ
        @inbounds for i in 1:length(basis.coeffs)
            !isassigned(basis.coeffs, i) && continue
            coeffs[i] = Vector{C}(undef, length(basis.coeffs[i]))
            for j in 1:length(basis.coeffs[i])
                coeffs[i][j] = basis.coeffs[i][j]
            end
        end
    else    # For Z and Q
        @inbounds for i in 1:length(basis.coeffs)
            !isassigned(basis.coeffs, i) && continue
            coeffs[i] = Vector{C}(undef, length(basis.coeffs[i]))
            for j in 1:length(basis.coeffs[i])
                # We cannot just use copy, since we mutate BigInts
                coeffs[i][j] = deepcopy(basis.coeffs[i][j])
            end
        end
    end

    basis_deep_copy_with_new_coeffs(basis, coeffs)
end

function basis_resize_if_needed!(basis::Basis{T}, to_add::Int) where {T}
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

# Normalize each element of the basis to have leading coefficient equal to 1
@timeit function basis_normalize!(
    basis::Basis{C},
    arithmetic::AbstractArithmeticZp{A, C}
) where {A <: Union{CoeffZp, CompositeCoeffZp}, C <: Union{CoeffZp, CompositeCoeffZp}}
    @log level = -5 "Normalizing polynomials in the basis"
    cfs = basis.coeffs
    @inbounds for i in 1:(basis.nfilled)
        !isassigned(cfs, i) && continue   # TODO: this is kind of bad
        mul = inv_mod_p(A(cfs[i][1]), arithmetic)
        cfs[i][1] = one(C)
        for j in 2:length(cfs[i])
            cfs[i][j] = mod_p(A(cfs[i][j]) * A(mul), arithmetic) % C
        end
        @invariant isone(cfs[i][1])
    end
    basis
end

# Normalize each element of the basis by dividing it by its leading coefficient
function basis_normalize!(
    basis::Basis{C},
    arithmetic::AbstractArithmeticQQ
) where {C <: CoeffQQ}
    @log level = -5 "Normalizing polynomials in the basis"
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
@timeit function pairset_update!(
    pairset::Pairset,
    basis::Basis{C},
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {C <: Coeff, M <: Monom}
    pr = monom_entrytype(M)
    pl, bl = pairset.load, idx
    ps = pairset.pairs
    lcms = pairset.lcms

    new_lead = basis.monoms[idx][1]

    # generate a pair for each pair
    @inbounds for i in 1:(bl - 1)
        newidx = pl + i
        if !basis.isredundant[i] &&
           !monom_is_gcd_const(ht.monoms[basis.monoms[i][1]], ht.monoms[new_lead])
            lcms[i] = get_lcm(basis.monoms[i][1], new_lead, ht, update_ht)
            deg = update_ht.hashdata[lcms[i]].deg
            ps[newidx] = CriticalPair(Int32(i), Int32(idx), lcms[i], pr(deg))
        else
            # lcm == 0 will mark redundancy of an S-pair
            lcms[i] = MonomId(0)
            # ps[newidx] = CriticalPair(i, idx, MonomId(0), pr(deg))
            ps[newidx] = CriticalPair{pr}(Int32(i), Int32(idx), MonomId(0), typemax(pr))
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
        if ps[i].deg > m && monom_is_divisible(ps[i].lcm, new_lead, ht)
            # mark an existing pair redundant
            ps[i] = CriticalPair{pr}(ps[i].poly1, ps[i].poly2, MonomId(0), ps[i].deg)
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
    insert_lcms_in_basis_hashtable!(pairset, pl, ht, update_ht, basis, lcms, j, pc + 1)

    # mark redundant polynomials in basis
    nonred = basis.nonredundant
    lml = basis.nnonredundant
    @inbounds for i in 1:lml
        if !basis.isredundant[nonred[i]]
            if monom_is_divisible(basis.monoms[nonred[i]][1], new_lead, ht)
                basis.isredundant[nonred[i]] = true
            end
        end
    end

    nothing
end

# Updates information about redundant generators in the basis
@timeit function basis_update!(basis::Basis, ht::MonomialHashtable{M}) where {M <: Monom}
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
function basis_is_new_polynomial_redundant!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {M}
    pt = monom_entrytype(M)
    resize_hashtable_if_needed!(update_ht, 0)

    # lead of new polynomial
    lead_new = basis.monoms[idx][1]
    ps = pairset.pairs

    @inbounds for i in (idx + 1):(basis.nfilled)
        i == idx && continue
        basis.isredundant[i] && continue

        # lead of new polynomial at index i > idx
        lead_i = basis.monoms[i][1]

        if monom_is_divisible(lead_new, lead_i, ht)
            # add new S-pair corresponding to Spoly(i, idx)
            lcm_new = get_lcm(lead_i, lead_new, ht, ht)
            psidx = pairset.load + 1
            ps[psidx] = CriticalPair{pt}(
                Int32(i),
                Int32(idx),
                lcm_new,
                pt(ht.hashdata[lcm_new].deg)
            )

            # mark redundant
            basis.isredundant[idx] = true
            pairset.load += 1

            return true
        end
    end

    return false
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
        basis.monoms[i] = Vector{MonomId}(undef, nterms)
        poly = basis.monoms[i]
        @inbounds for j in 1:nterms
            poly[j] = insert_in_hashtable!(ht, exponents[i][j])
        end
    end

    basis.nfilled = ngens
end

function basis_sweep_redundant!(basis::Basis, hashtable)
    # here -- assert that basis is in fact a Groebner basis.
    # NOTE: maybe sort generators for more effective sweeping?
    @inbounds for i in 1:(basis.nprocessed)
        for j in (i + 1):(basis.nprocessed)
            basis.isredundant[i] && continue
            basis.isredundant[j] && continue
            lead_i = basis.monoms[i][1]
            lead_j = basis.monoms[j][1]
            if monom_is_divisible(lead_i, lead_j, hashtable)
                basis.isredundant[i] = true
            elseif monom_is_divisible(lead_j, lead_i, hashtable)
                basis.isredundant[j] = true
            end
        end
    end
    nothing
end

# Remove redundant elements from the basis by moving all non-redundant up front
function basis_mark_redundant_elements!(basis::Basis)
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

@timeit function basis_standardize!(
    ring,
    basis::Basis,
    ht::MonomialHashtable,
    ord,
    arithmetic
)
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
    basis_normalize!(basis, arithmetic)
end

# Returns the monomials of the polynomials in the basis
function basis_get_monoms_by_identifiers(
    basis::Basis,
    ht::MonomialHashtable{M}
) where {M <: Monom}
    monoms = Vector{Vector{M}}(undef, basis.nnonredundant)
    @inbounds for i in 1:(basis.nnonredundant)
        idx = basis.nonredundant[i]
        poly = basis.monoms[idx]
        monoms[i] = Vector{M}(undef, length(poly))
        for j in 1:length(poly)
            monoms[i][j] = ht.monoms[poly[j]]
        end
    end
    monoms
end

# Returns the monomials and the coefficients of polynomials in the basis
@timeit function basis_export_data(
    basis::Basis{C},
    ht::MonomialHashtable{M}
) where {M <: Monom, C <: Coeff}
    exps = basis_get_monoms_by_identifiers(basis, ht)
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
function insert_lcms_in_basis_hashtable!(
    pairset::Pairset,
    off::Int,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    basis::Basis,
    plcm::Vector{MonomId},
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

        if monom_is_gcd_const(
            monoms[ps[off + l].poly1][1],
            monoms[ps[off + 1].poly2][1],
            ht
        )
            l += 1
            continue
        end

        ps[m] = ps[off + l]

        h = update_ht.hashdata[plcm[l]].hash
        ht.monoms[ht.load + 1] = monom_copy(update_ht.monoms[plcm[l]])
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

            ps[m] = CriticalPair{typeof(ps[m].deg)}(ps[m].poly1, ps[m].poly2, hm, ps[m].deg)
            m += 1
            l += 1
            @goto Letsgo
        end

        @invariant !ht.frozen

        ht.hashtable[k] = pos = ht.load + 1

        uhd = update_ht.hashdata
        ll = plcm[l]
        ht.hashdata[ht.load + 1] = Hashvalue(0, h, uhd[ll].divmask, uhd[ll].deg)

        ht.load += 1
        ps[m] = CriticalPair{typeof(ps[m].deg)}(
            ps[m].poly1,
            ps[m].poly2,
            MonomId(pos),
            ps[m].deg
        )
        m += 1
        l += 1
    end

    pairset.load = m - 1
end
