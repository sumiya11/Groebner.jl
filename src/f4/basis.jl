# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve
#   https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+
#   https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# Pairset and Basis

# Pairset is a list of critical pairs (CriticalPair).

# Basis is a structure that stores a list of polynomials. Each polynomial is
# represented with a sorted vector of monomials and a vector of coefficients.
# Monomials and coefficients are stored in the basis separately. Each monomial
# is represented with an integer -- a unique identifier that indexes a bucket in
# the hashtable.

###
# Pairset

const CRITICAL_PAIR_REDUNDANT = MonomId(0)

# A pair of polynomials
struct CriticalPair
    # First polynomial given by its index in the basis array
    poly1::Int32
    # Second polynomial given by its index in the basis array
    poly2::Int32
    # Index of lcm(lead(poly1), lead(poly2)) in the hashtable
    lcm::MonomId
end

# Stores S-Pairs and some additional info.
mutable struct Pairset{ExponentType <: Integer}
    pairs::Vector{CriticalPair}
    degrees::Vector{ExponentType}   # Buffer for lcms' degrees
    lcms::Vector{MonomId}           # Buffer for lcms
    load::Int
    scratch1::Vector{Int}            # Scratch spaces for faster sorting
    scratch2::Vector{CriticalPair}
    scratch3::Vector{ExponentType}
end

function pairset_initialize(::Type{ExponentType}; initial_size=2^6) where {ExponentType}
    Pairset(
        Vector{CriticalPair}(undef, initial_size),
        Vector{ExponentType}(undef, initial_size),
        Vector{MonomId}(),
        0,
        Vector{Int}(),
        Vector{CriticalPair}(),
        Vector{ExponentType}()
    )
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
    resize!(ps.degrees, newsize)
    nothing
end

function pairset_resize_lcms_if_needed!(ps::Pairset, nfilled::Int)
    if length(ps.lcms) < nfilled + 1
        # NOTE: Resizing by a small factor is questionable.
        resize!(ps.lcms, floor(Int, nfilled * 1.1) + 1)
    end
    nothing
end

function pairset_find_smallest_degree_pair(ps::Pairset)
    @invariant ps.load > 0 && length(ps.pairs) > 0 && length(ps.degrees) > 0
    degs = ps.degrees
    @inbounds pair_idx, pair_min_deg = 1, degs[1]
    @inbounds for i in 1:(ps.load)
        if degs[i] <= pair_min_deg
            pair_min_deg = degs[i]
            pair_idx = i
        end
    end
    pair_idx, pair_min_deg
end

function pairset_partition_by_degree!(ps::Pairset)
    @invariant ps.load > 0
    _, pair_min_deg = pairset_find_smallest_degree_pair(ps)

    pairs = ps.pairs
    degs = ps.degrees
    i, j = 0, ps.load + 1
    @inbounds while true
        i += 1
        j -= 1
        while i <= ps.load && degs[i] == pair_min_deg
            i += 1
        end
        while j > 1 && degs[j] > pair_min_deg
            j -= 1
        end
        i >= j && break
        pairs[i], pairs[j] = pairs[j], pairs[i]
        degs[i], degs[j] = degs[j], degs[i]
    end

    i - 1
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

    changematrix::Vector{Dict{Int, Dict{MonomId, C}}}
end

# Initialize basis with coefficient of type T.
function basis_initialize(ring::PolyRing, sz::Int, ::Type{T}) where {T <: Coeff}
    Basis(
        Vector{Vector{MonomId}}(undef, sz),
        Vector{Vector{T}}(undef, sz),
        sz,
        0,
        0,
        zeros(Bool, sz),
        Vector{Int}(undef, sz),
        Vector{DivisionMask}(undef, sz),
        0,
        Vector{Int}(undef, 0),
        Vector{Dict{Int, Dict{MonomId, T}}}(undef, 0)
    )
end

# Initialize basis with the given (already hashed) monomials and coefficients.
function basis_initialize(
    ring::PolyRing,
    hashedexps::Vector{Vector{MonomId}},
    coeffs::Vector{Vector{T}}
) where {T <: Coeff}
    sz = length(hashedexps)
    Basis(
        hashedexps,
        coeffs,
        sz,
        0,
        sz,
        zeros(Bool, sz),
        Vector{Int}(undef, sz),
        Vector{DivisionMask}(undef, sz),
        0,
        Vector{Int}(undef, 0),
        Vector{Dict{Int, Dict{MonomId, T}}}(undef, 0)
    )
end

# Same as basis_initialize, but uses an existing hashtable
function basis_initialize_using_existing_hashtable(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable;
) where {M, C <: Coeff}
    basis = basis_initialize(ring, length(monoms), C)
    basis_fill_data!(basis, present_ht, monoms, coeffs)
    basis
end

###
# Change matrix

function basis_changematrix_initialize!(
    basis::Basis{C},
    hashtable::MonomialHashtable{M, Ord}
) where {C <: Coeff, M <: Monom, Ord}
    resize!(basis.changematrix, basis.nfilled)
    id_of_1 = hashtable_insert!(hashtable, monom_construct_const(M, hashtable.nvars))
    for i in 1:(basis.nfilled)
        basis.changematrix[i] = Dict{Int, Dict{MonomId, C}}()
        basis.changematrix[i][i] = Dict{MonomId, C}(id_of_1 => one(C))
    end
    nothing
end

function basis_changematrix_mul!(
    basis::Basis,
    idx::Int,
    cf::C,
    arithmetic::AbstractArithmetic
) where {C <: Coeff}
    row = basis.changematrix[idx]
    new_row = empty(row)
    for (pos, poly) in row
        new_poly = empty(poly)
        for (monom_id, monom_cf) in poly
            new_poly[monom_id] = mod_p(monom_cf * cf, arithmetic)
        end
        new_row[pos] = new_poly
    end
    basis.changematrix[idx] = new_row
    nothing
end

function basis_changematrix_addmul!(
    basis::Basis,
    ht::MonomialHashtable,
    symbol_ht,
    idx::Int,
    poly_idx,
    poly_mult,
    cf::C,
    arithmetic::AbstractArithmeticZp{AccumType, CoeffType}
) where {C <: Coeff, AccumType, CoeffType}
    row = basis.changematrix[idx]
    if true
        @inbounds for (ref_poly_idx, ref_quo) in basis.changematrix[poly_idx]
            if !haskey(row, ref_poly_idx)
                row[ref_poly_idx] = Dict{MonomId, C}()
            end
            poly = row[ref_poly_idx]
            hashtable_resize_if_needed!(ht, length(ref_quo))
            for (ref_mult, ref_cf) in ref_quo
                ref_monom = ht.monoms[ref_mult]
                quo_monom = ht.monoms[poly_mult]
                new_monom = monom_copy(ref_monom)
                new_monom = monom_product!(new_monom, ref_monom, quo_monom)
                new_monom_id = hashtable_insert!(ht, new_monom)
                if !haskey(poly, new_monom_id)
                    poly[new_monom_id] = zero(CoeffType)
                end
                poly[new_monom_id] =
                    mod_p(poly[new_monom_id] + ref_cf * AccumType(cf), arithmetic)
            end
        end
    else
        if !haskey(row, poly_idx)
            row[poly_idx] = Dict{MonomId, C}()
        end
        poly = row[poly_idx]
        if !haskey(poly, poly_mult)
            poly[poly_mult] = zero(CoeffType)
        end
        poly[poly_mult] = mod_p(poly[poly_mult] + AccumType(cf), arithmetic)
    end
    nothing
end

function basis_changematrix_deep_copy_with_new_type(
    changematrix,
    new_sparse_row_coeffs::Vector{Vector{C}}
) where {C}
    new_changematrix = Vector{Dict{Int, Dict{MonomId, C}}}(undef, length(changematrix))
    for i in 1:length(changematrix)
        !isassigned(changematrix, i) && continue
        obj = changematrix[i]
        new_obj = Dict{Int, Dict{MonomId, C}}()
        for (k, v) in obj
            new_obj[k] = Dict{MonomId, C}()
            for (k2, v2) in v
                new_obj[k][k2] = C(v2)
            end
        end
        new_changematrix[i] = new_obj
    end
    new_changematrix
end

function basis_changematrix_shallow_copy_with_new_type(changematrix, new_sparse_row_coeffs)
    basis_changematrix_deep_copy_with_new_type(changematrix, new_sparse_row_coeffs)
end

function basis_changematrix_export(
    basis::Basis{C},
    ht::MonomialHashtable{M},
    npolys
) where {C <: Coeff, M <: Monom}
    @log :all basis.changematrix
    matrix_monoms = Vector{Vector{Vector{M}}}(undef, basis.nnonredundant)
    matrix_coeffs = Vector{Vector{Vector{C}}}(undef, basis.nnonredundant)
    @inbounds for i in 1:(basis.nnonredundant)
        matrix_monoms[i] = Vector{Vector{M}}(undef, npolys)
        matrix_coeffs[i] = Vector{Vector{M}}(undef, npolys)
        idx = basis.nonredundant[i]
        row = basis.changematrix[idx]
        for poly_idx in 1:npolys
            matrix_monoms[i][poly_idx] = Vector{M}()
            matrix_coeffs[i][poly_idx] = Vector{C}()
            !haskey(row, poly_idx) && continue
            for (monom_id, monom_cf) in row[poly_idx]
                push!(matrix_monoms[i][poly_idx], ht.monoms[monom_id])
                push!(matrix_coeffs[i][poly_idx], monom_cf)
            end
        end
        sort_input_terms_to_change_ordering!(matrix_monoms[i], matrix_coeffs[i], ht.ord)
    end
    matrix_monoms, matrix_coeffs
end

###
# Basis utils

function basis_shallow_copy_with_new_coeffs(
    basis::Basis{C},
    new_sparse_row_coeffs::Vector{Vector{T}}
) where {C <: Coeff, T <: Coeff}
    Basis(
        basis.monoms,
        new_sparse_row_coeffs,
        basis.size,
        basis.nprocessed,
        basis.nfilled,
        basis.isredundant,
        basis.nonredundant,
        basis.divmasks,
        basis.nnonredundant,
        basis.sugar_cubes,
        basis_changematrix_shallow_copy_with_new_type(
            basis.changematrix,
            new_sparse_row_coeffs
        )
    )
end

function basis_deep_copy_with_new_coeffs(
    basis::Basis{C},
    new_sparse_row_coeffs::Vector{Vector{T}}
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
        new_sparse_row_coeffs,
        basis.size,
        basis.nprocessed,
        basis.nfilled,
        copy(basis.isredundant),
        copy(basis.nonredundant),
        copy(basis.divmasks),
        basis.nnonredundant,
        copy(basis.sugar_cubes),
        basis_changematrix_deep_copy_with_new_type(
            basis.changematrix,
            new_sparse_row_coeffs
        )
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
        # resize!(basis.sugar_cubes, basis.size)
    end
    @invariant basis.size >= basis.nprocessed + to_add
    nothing
end

# Normalize each element of the basis to have leading coefficient equal to 1
function basis_make_monic!(
    basis::Basis{C},
    arithmetic::AbstractArithmeticZp{A, C},
    changematrix::Bool
) where {A <: Union{CoeffZp, CompositeCoeffZp}, C <: Union{CoeffZp, CompositeCoeffZp}}
    @log :debug "Normalizing polynomials in the basis"
    cfs = basis.coeffs
    @inbounds for i in 1:(basis.nfilled)
        !isassigned(cfs, i) && continue   # TODO: this is kind of bad
        mul = inv_mod_p(A(cfs[i][1]), arithmetic)
        cfs[i][1] = one(C)
        for j in 2:length(cfs[i])
            cfs[i][j] = mod_p(A(cfs[i][j]) * A(mul), arithmetic) % C
        end
        @invariant isone(cfs[i][1])
        if changematrix
            basis_changematrix_mul!(basis, i, A(mul), arithmetic)
        end
    end
    basis
end

# Normalize each element of the basis by dividing it by its leading coefficient
function basis_make_monic!(
    basis::Basis{C},
    arithmetic::AbstractArithmeticQQ,
    changematrix::Bool
) where {C <: CoeffQQ}
    @log :debug "Normalizing polynomials in the basis"
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
@timeit function pairset_update!(
    pairset::Pairset{D},
    basis::Basis{C},
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {D, C <: Coeff, M <: Monom}
    pl, bl = pairset.load, idx
    ps = pairset.pairs
    lcms = pairset.lcms
    degs = pairset.degrees

    @inbounds new_lead = basis.monoms[idx][1]

    # Generate new pairs.
    @inbounds for i in 1:(bl - 1)
        newidx = pl + i
        if !basis.isredundant[i] &&
           !monom_is_gcd_const(ht.monoms[basis.monoms[i][1]], ht.monoms[new_lead])
            lcms[i] = hashtable_get_lcm!(basis.monoms[i][1], new_lead, ht, update_ht)
            degs[newidx] = monom_totaldeg(update_ht.monoms[lcms[i]])
            ps[newidx] = CriticalPair(Int32(i), Int32(idx), lcms[i])
        else
            lcms[i] = CRITICAL_PAIR_REDUNDANT
            degs[newidx] = typemax(D)
            ps[newidx] = CriticalPair(Int32(i), Int32(idx), CRITICAL_PAIR_REDUNDANT)
        end
    end

    # Traverse existing pairs...
    @inbounds for i in 1:pl
        if ps[i].lcm == CRITICAL_PAIR_REDUNDANT
            continue
        end

        j = ps[i].poly1
        l = ps[i].poly2
        m = max(degs[pl + l], degs[pl + j])

        # ...if an existing pair is divisible by the lead of the new poly and
        # has a greater degree than newly generated critical pair, then mark the
        # existing pair redundant
        if degs[i] > m && hashtable_monom_is_divisible(ps[i].lcm, new_lead, ht)
            ps[i] = CriticalPair(ps[i].poly1, ps[i].poly2, CRITICAL_PAIR_REDUNDANT)
        end
    end

    # Traverse new pairs to move non-redundant ones first.
    j = 1
    @inbounds for i in 1:(bl - 1)
        if !basis.isredundant[i]
            ps[pl + j] = ps[pl + i]
            degs[pl + j] = degs[pl + i]
            j += 1
        end
    end

    sort_pairset_by_degree!(pairset, pl + 1, j - 2)

    @inbounds for i in 1:(j - 1)
        lcms[i] = ps[pl + i].lcm
    end
    @inbounds lcms[j] = CRITICAL_PAIR_REDUNDANT
    pc = j
    pc -= 1

    # Mark redundancy of some pairs based on their lcms.
    @inbounds for j in 1:pc
        if !(lcms[j] == CRITICAL_PAIR_REDUNDANT)
            hashtable_check_monomial_division_in_update(lcms, j + 1, pc, lcms[j], update_ht)
        end
    end

    # Remove redundant pairs from the pairset.
    j = 1
    @inbounds for i in 1:(pairset.load)
        (ps[i].lcm == CRITICAL_PAIR_REDUNDANT) && continue
        ps[j] = ps[i]
        degs[j] = degs[i]
        j += 1
    end

    hashtable_resize_if_needed!(ht, pc)

    # Add new lcm monomials to the basis hashtable 
    # (including index j and not including index pc).
    insert_lcms_in_basis_hashtable!(pairset, pl, ht, update_ht, basis, lcms, j, pc + 1)

    # Mark redundant polynomials in basis.
    nonred = basis.nonredundant
    lml = basis.nnonredundant
    @inbounds for i in 1:lml
        if !basis.isredundant[nonred[i]]
            if hashtable_monom_is_divisible(basis.monoms[nonred[i]][1], new_lead, ht)
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

# Checks if the element of basis at position idx is redundant.
function basis_is_new_polynomial_redundant!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {M <: Monom}
    hashtable_resize_if_needed!(update_ht, 0)

    @inbounds lead_new = basis.monoms[idx][1]
    ps = pairset.pairs
    degs = pairset.degrees

    @inbounds for i in (idx + 1):(basis.nfilled)
        basis.isredundant[i] && continue

        lead_i = basis.monoms[i][1]

        !hashtable_monom_is_divisible(lead_new, lead_i, ht) && continue

        # Add a new critical pair corresponding to Spoly(i, idx).
        pairset_resize_if_needed!(pairset, 1)
        lcm_new = hashtable_get_lcm!(lead_i, lead_new, ht, ht)
        ps[pairset.load + 1] = CriticalPair(Int32(i), Int32(idx), lcm_new)
        degs[pairset.load + 1] = monom_totaldeg(ht.monoms[lcm_new])
        pairset.load += 1

        # Mark the polynomial as redundant.
        basis.isredundant[idx] = true
        return true
    end

    false
end

# given input exponent and coefficient vectors hashes exponents into `ht`
# and then constructs hashed polynomials for `basis`
function basis_fill_data!(
    basis::Basis,
    ht::MonomialHashtable{M},
    exponents::Vector{Vector{M}},
    coeffs::Vector{Vector{T}}
) where {M, T}
    ngens = length(exponents)
    @inbounds for i in 1:ngens
        hashtable_resize_if_needed!(ht, length(exponents[i]))

        nterms = length(coeffs[i])
        basis.coeffs[i] = coeffs[i]
        basis.monoms[i] = Vector{MonomId}(undef, nterms)
        poly = basis.monoms[i]
        @inbounds for j in 1:nterms
            poly[j] = hashtable_insert!(ht, exponents[i][j])
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
            if hashtable_monom_is_divisible(lead_i, lead_j, hashtable)
                basis.isredundant[i] = true
            elseif hashtable_monom_is_divisible(lead_j, lead_i, hashtable)
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
    @invariant basis.nprocessed == basis.nfilled
    basis
end

@timeit function basis_standardize!(
    ring,
    basis::Basis,
    ht::MonomialHashtable,
    ord,
    arithmetic,
    changematrix::Bool
)
    @inbounds for i in 1:(basis.nnonredundant)
        idx = basis.nonredundant[i]
        basis.nonredundant[i] = i
        basis.isredundant[i] = false
        basis.coeffs[i] = basis.coeffs[idx]
        basis.monoms[i] = basis.monoms[idx]
        if changematrix
            basis.changematrix[i] = basis.changematrix[idx]
        end
    end
    basis.size = basis.nprocessed = basis.nfilled = basis.nnonredundant
    resize!(basis.coeffs, basis.nprocessed)
    resize!(basis.monoms, basis.nprocessed)
    resize!(basis.divmasks, basis.nprocessed)
    resize!(basis.nonredundant, basis.nprocessed)
    resize!(basis.isredundant, basis.nprocessed)
    resize!(basis.changematrix, basis.nprocessed)
    # resize!(basis.sugar_cubes, basis.nprocessed)
    sort_polys_by_lead_increasing!(basis, ht, changematrix, ord=ord)
    basis_make_monic!(basis, arithmetic, changematrix)
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
    degs = pairset.degrees

    mod = MonomHash(ht.size - 1)
    @invariant ispow2(mod + 1)

    m = ifirst
    l = 1
    @label Letsgo
    @inbounds while l < ilast
        if plcm[l] == CRITICAL_PAIR_REDUNDANT
            l += 1
            continue
        end

        if hashtable_monom_is_gcd_const(
            monoms[ps[off + l].poly1][1],
            monoms[ps[off + 1].poly2][1],
            ht
        )
            l += 1
            continue
        end

        ps[m] = ps[off + l]
        degs[m] = degs[off + l]

        h = update_ht.hashdata[plcm[l]].hash
        ht.monoms[ht.load + 1] = monom_copy(update_ht.monoms[plcm[l]])
        n = ht.monoms[ht.load + 1]

        k = h
        i = MonomHash(0)
        @inbounds while i <= ht.size
            k = hashtable_next_lookup_index(h, i, mod)
            hm = ht.hashtable[k]

            # if free
            iszero(hm) && break

            if hashtable_is_hash_collision(ht, hm, n, h)
                i += MonomHash(1)
                continue
            end

            ps[m] = CriticalPair(ps[m].poly1, ps[m].poly2, hm)
            m += 1
            l += 1
            @goto Letsgo
        end

        @invariant !ht.frozen

        ht.hashtable[k] = pos = ht.load + 1

        uhd = update_ht.hashdata
        ll = plcm[l]
        ht.hashdata[ht.load + 1] = Hashvalue(0, h, uhd[ll].divmask)

        ht.load += 1
        ps[m] = CriticalPair(ps[m].poly1, ps[m].poly2, MonomId(pos))
        m += 1
        l += 1
    end

    pairset.load = m - 1
end
