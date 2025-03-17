# This file is a part of Groebner.jl. License is GNU GPL v2.

# Parts of this file were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING

###
# Pairset

const CRITICAL_PAIR_REDUNDANT = MonomId(0)

struct CriticalPair
    poly1::Int32
    poly2::Int32
    # Index of lcm(lm(poly1), lm(poly2)) in the hashtable
    lcm::MonomId

    function CriticalPair(poly1, poly2, lcm)
        @invariant poly1 < poly2
        new(poly1, poly2, lcm)
    end
end

mutable struct Pairset{ExponentType <: Integer}
    pairs::Vector{CriticalPair}
    degrees::Vector{ExponentType}
    lcms::Vector{MonomId}
    load::Int
    scratch1::Vector{Int}           # Scratch spaces for faster sorting
    scratch2::Vector{CriticalPair}
    scratch3::Vector{ExponentType}
    scratch4::Vector{Int}
end

function pairset_initialize(::Type{ExponentType}; initial_size=2^6) where {ExponentType}
    Pairset(
        Vector{CriticalPair}(undef, initial_size),
        Vector{ExponentType}(undef, initial_size),
        Vector{MonomId}(),
        0,
        Vector{Int}(),
        Vector{CriticalPair}(),
        Vector{ExponentType}(),
        Vector{Int}()
    )
end

Base.isempty(ps::Pairset) = ps.load == 0

function pairset_resize_if_needed!(ps::Pairset, to_add::Int)
    newsize = length(ps.pairs)
    while ps.load + to_add > newsize
        newsize = max(2 * newsize, ps.load + to_add)
    end
    resize!(ps.pairs, newsize)
    resize!(ps.degrees, newsize)
    nothing
end

function pairset_resize_lcms_if_needed!(ps::Pairset, n_filled::Int)
    if length(ps.lcms) < n_filled + 1
        resize!(ps.lcms, floor(Int, n_filled * 1.1) + 1)
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

# Returns N, the number of critical pairs of the smallest degree. Sorts the
# critical pairs so that the first N pairs in the pairset are the smallest with
# respect to degree.
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

# Basis is a list of polynomials. A polynomial is represented by a sorted vector
# of monomials and a sorted vector of coefficients. Monomials and coefficients
# are stored in the basis separately. A monomial is represented by an integer, a
# unique identifier provided by the hashtable.
mutable struct Basis{C <: Coeff}
    # Monomial identifiers are integers greater than zero. 
    # Zero and negative values are reserved.
    monoms::Vector{Vector{MonomId}}
    coeffs::Vector{Vector{C}}

    # filled >= processed at any time.
    n_filled::Int
    n_processed::Int

    n_nonredundant::Int
    is_redundant::Vector{Bool}
    nonredundant_indices::Vector{Int}

    divmasks::Vector{DivisionMask}

    # usually empty
    changematrix::Vector{Dict{Int, Dict{MonomId, C}}}
end

function basis_initialize(ring::PolyRing, sz::Int, ::Type{C}) where {C <: Coeff}
    Basis(
        Vector{Vector{MonomId}}(undef, sz),
        Vector{Vector{C}}(undef, sz),
        0,
        0,
        0,
        zeros(Bool, sz),
        Vector{Int}(undef, sz),
        Vector{DivisionMask}(undef, sz),
        Vector{Dict{Int, Dict{MonomId, C}}}(undef, 0)
    )
end

# Same as basis_initialize, but uses an existing hashtable.
function basis_initialize_using_existing_hashtable(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    present_ht::MonomialHashtable;
) where {M <: Monom, C <: Coeff}
    basis = basis_initialize(ring, length(monoms), C)
    basis_fill_data!(basis, present_ht, monoms, coeffs)
    basis
end

function basis_well_formed(ring::PolyRing, basis::Basis, hashtable::MonomialHashtable)
    (isempty(basis.monoms) || isempty(basis.coeffs)) && error("Basis cannot be empty")
    !(basis.n_filled >= basis.n_processed) && error("Basis cannot be empty")
    !is_sorted_by_lead_increasing(basis, hashtable) && error("Basis elements must be sorted")
    for i in basis.n_filled
        isempty(basis.monoms[i]) && error("Zero polynomials are not allowed.")
        !(length(basis.monoms[i]) == length(basis.coeffs[i])) && error("Bad polynomial")
        !allunique(basis.monoms[i]) && error("Bad polynomial")
        !issorted(
            basis.monoms[i],
            lt=(j, k) -> monom_isless(hashtable.monoms[k], hashtable.monoms[j], ring.ord)
        ) && error("Bad polynomial")
	iszero(basis.coeffs[i][1]) && error("Leading coefficient is zero")
	for j in 1:length(basis.coeffs[i])
            (ring.ground == :zp) &&
	    !(-1 < basis.coeffs[i][j] < ring.characteristic) &&
                error("Coefficients must be normalized")
            (basis.monoms[i][j] > hashtable.load) && error("Bad monomial")
        end
    end
    true
end

###
# Change matrix

function basis_changematrix_initialize!(
    basis::Basis{C},
    hashtable::MonomialHashtable{M, Ord}
) where {C <: Coeff, M <: Monom, Ord}
    resize!(basis.changematrix, basis.n_filled)
    id_of_1 = hashtable_insert!(hashtable, monom_construct_const(M, hashtable.nvars))
    for i in 1:(basis.n_filled)
        basis.changematrix[i] = Dict{Int, Dict{MonomId, C}}()
        basis.changematrix[i][i] = Dict{MonomId, C}(id_of_1 => one(C))
    end
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
            poly[new_monom_id] = mod_p(poly[new_monom_id] + ref_cf * AccumType(cf), arithmetic)
        end
    end
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
    matrix_monoms = Vector{Vector{Vector{M}}}(undef, basis.n_nonredundant)
    matrix_coeffs = Vector{Vector{Vector{C}}}(undef, basis.n_nonredundant)
    @inbounds for i in 1:(basis.n_nonredundant)
        matrix_monoms[i] = Vector{Vector{M}}(undef, npolys)
        matrix_coeffs[i] = Vector{Vector{M}}(undef, npolys)
        idx = basis.nonredundant_indices[i]
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
        basis.n_filled,
        basis.n_processed,
        basis.n_nonredundant,
        basis.is_redundant,
        basis.nonredundant_indices,
        basis.divmasks,
        basis_changematrix_shallow_copy_with_new_type(basis.changematrix, new_sparse_row_coeffs)
    )
end

function basis_deep_copy_with_new_coeffs(
    basis::Basis{C},
    new_sparse_row_coeffs::Vector{Vector{T}}
) where {C <: Coeff, T <: Coeff}
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
        basis.n_filled,
        basis.n_processed,
        basis.n_nonredundant,
        copy(basis.is_redundant),
        copy(basis.nonredundant_indices),
        copy(basis.divmasks),
        basis_changematrix_deep_copy_with_new_type(basis.changematrix, new_sparse_row_coeffs)
    )
end

function basis_deepcopy(basis::Basis{C}) where {C <: Coeff}
    coeffs = Vector{Vector{C}}(undef, length(basis.coeffs))

    if isbitstype(C)  # For Z/pZ
        @inbounds for i in 1:length(basis.coeffs)
            !isassigned(basis.coeffs, i) && continue
            coeffs[i] = Vector{C}(undef, length(basis.coeffs[i]))
            for j in 1:length(basis.coeffs[i])
                coeffs[i][j] = basis.coeffs[i][j]
            end
        end
    else  # For Z and Q
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
    size = length(basis.monoms)
    while basis.n_processed + to_add >= size
        size = max(size * 2, basis.n_processed + to_add)
        resize!(basis.monoms, size)
        resize!(basis.coeffs, size)
        resize!(basis.is_redundant, size)
        @inbounds basis.is_redundant[(basis.n_processed + 1):end] .= false
        resize!(basis.nonredundant_indices, size)
        resize!(basis.divmasks, size)
    end
    @invariant size >= basis.n_processed + to_add
    nothing
end

function basis_make_monic!(
    basis::Basis{C},
    arithmetic::AbstractArithmeticZp{A, C},
    changematrix::Bool
) where {A <: Union{CoeffZp, CompositeCoeffZp}, C <: Union{CoeffZp, CompositeCoeffZp}}
    cfs = basis.coeffs
    @inbounds for i in 1:(basis.n_filled)
        !isassigned(cfs, i) && continue
        isone(cfs[i][1]) && continue
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

function basis_make_monic!(
    basis::Basis{C},
    arithmetic::AbstractArithmetic,
    changematrix::Bool
) where {C <: Union{CoeffQQ, CoeffGeneric}}
    cfs = basis.coeffs
    @inbounds for i in 1:(basis.n_filled)
        !isassigned(cfs, i) && continue
        isone(cfs[i][1]) && continue
        mul = inv(cfs[i][1])
        for j in 2:length(cfs[i])
            cfs[i][j] *= mul
        end
        cfs[i][1] = one(cfs[i][1])
    end
    basis
end

# Follows paragraph 4.4. from 
# Gebauer, Möller, On an Installation of Buchberger's Algorithm
function pairset_update!(
    pairset::Pairset{D},
    basis::Basis{C},
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {D, C <: Coeff, M <: Monom}
    pairset_resize_lcms_if_needed!(pairset, basis.n_filled)
    pairset_resize_if_needed!(pairset, basis.n_filled)

    pairs = pairset.pairs
    lcms = pairset.lcms
    degs = pairset.degrees

    @inbounds lead_idx = basis.monoms[idx][1]

    # Generate new pairs.
    hashtable_resize_if_needed!(update_ht, idx)
    @inbounds for i in 1:(idx - 1)
        newidx = pairset.load + i
        if !basis.is_redundant[i] &&
           !monom_is_gcd_const(ht.monoms[basis.monoms[i][1]], ht.monoms[lead_idx])
            lcms[i] = hashtable_get_lcm!(basis.monoms[i][1], lead_idx, ht, update_ht)
            degs[newidx] = monom_totaldeg(update_ht.monoms[lcms[i]])
            pairs[newidx] = CriticalPair(Int32(i), Int32(idx), lcms[i])
        else
            lcms[i] = CRITICAL_PAIR_REDUNDANT
            degs[newidx] = typemax(D)
            pairs[newidx] = CriticalPair(Int32(i), Int32(idx), CRITICAL_PAIR_REDUNDANT)
        end
    end

    # Criterion B_k(j, l).
    # T(j, l) is divisible by T(k)
    # max(deg T(j, k), deg T(l, k)) < deg T(j, l)
    @inbounds for i in 1:(pairset.load)
        (pairs[i].lcm == CRITICAL_PAIR_REDUNDANT) && continue
        j = pairs[i].poly1
        l = pairs[i].poly2
        m = max(degs[pairset.load + l], degs[pairset.load + j])
        if degs[i] > m && hashtable_monom_is_divisible(pairs[i].lcm, lead_idx, ht)
            pairs[i] = CriticalPair(pairs[i].poly1, pairs[i].poly2, CRITICAL_PAIR_REDUNDANT)
        end
    end

    # Traverse new pairs to move non-redundant ones first.
    cnt = 1
    @inbounds for i in 1:(idx - 1)
        (lcms[i] == CRITICAL_PAIR_REDUNDANT) && continue
        # if !basis.is_redundant[i]
        pairs[pairset.load + cnt] = pairs[pairset.load + i]
        degs[pairset.load + cnt] = degs[pairset.load + i]
        cnt += 1
        # end
    end

    sort_pairset_by_degree!(pairset, pairset.load + 1, cnt - 2)

    @inbounds for i in 1:(cnt - 1)
        lcms[i] = pairs[pairset.load + i].lcm
    end
    @inbounds lcms[cnt] = CRITICAL_PAIR_REDUNDANT
    new_pairs = cnt - 1

    # Criterion M(i, j).
    # T(i, j) | T(k, j)
    @inbounds for i in 1:new_pairs
        (lcms[i] == CRITICAL_PAIR_REDUNDANT) && continue
        divmask_i = update_ht.divmasks[lcms[i]]
        lcm_i = update_ht.monoms[lcms[i]]
        k = i + 1
        for k in (i + 1):new_pairs
            (lcms[k] == CRITICAL_PAIR_REDUNDANT) && continue
            if update_ht.use_divmask &&
               !divmask_is_probably_divisible(update_ht.divmasks[lcms[k]], divmask_i)
                continue
            end
            ea = update_ht.monoms[lcms[k]]
            if monom_is_divisible(update_ht.monoms[lcms[k]], lcm_i)
                lcms[k] = CRITICAL_PAIR_REDUNDANT
            end
        end
    end

    # Remove redundant pairs from the pairset.
    cnt = 1
    @inbounds for i in 1:(pairset.load)
        (pairs[i].lcm == CRITICAL_PAIR_REDUNDANT) && continue
        pairs[cnt] = pairs[i]
        degs[cnt] = degs[i]
        cnt += 1
    end

    # Add new lcm monomials to the main hashtable.
    hashtable_resize_if_needed!(ht, new_pairs)
    @inbounds for i in 1:new_pairs
        (lcms[i] == CRITICAL_PAIR_REDUNDANT) && continue
        id = hashtable_insert!(ht, update_ht.monoms[lcms[i]])
        pair = pairs[pairset.load + i]
        pairs[cnt] = CriticalPair(pair.poly1, pair.poly2, id)
        degs[cnt] = degs[pairset.load + i]
        cnt += 1
    end
    pairset.load = cnt - 1

    # Mark redundant polynomials in the generating set.
    nonred = basis.nonredundant_indices
    @inbounds for i in 1:(basis.n_nonredundant)
        basis.is_redundant[nonred[i]] && continue
        if hashtable_monom_is_divisible(basis.monoms[nonred[i]][1], lead_idx, ht)
            basis.is_redundant[nonred[i]] = true
        end
    end
end

function basis_update!(basis::Basis, ht::MonomialHashtable{M}) where {M <: Monom}
    k = 1
    lead = basis.divmasks
    nonred = basis.nonredundant_indices
    @inbounds for i in 1:(basis.n_nonredundant)
        if !basis.is_redundant[nonred[i]]
            basis.divmasks[k] = lead[i]
            basis.nonredundant_indices[k] = nonred[i]
            k += 1
        end
    end
    basis.n_nonredundant = k - 1

    @inbounds for i in (basis.n_processed + 1):(basis.n_filled)
        if !basis.is_redundant[i]
            lead[k] = ht.divmasks[basis.monoms[i][1]]
            nonred[k] = i
            k += 1
        end
    end

    basis.n_nonredundant = k - 1
    basis.n_processed = basis.n_filled
end

function basis_mark_redundant_elements!(
    pairset::Pairset,
    basis::Basis,
    ht::MonomialHashtable{M},
    update_ht::MonomialHashtable{M},
    idx::Int
) where {M <: Monom}
    @inbounds lead_idx = basis.monoms[idx][1]
    ps = pairset.pairs
    degs = pairset.degrees
    @inbounds for i in (idx + 1):(basis.n_filled)
        basis.is_redundant[i] && continue
        lead_i = basis.monoms[i][1]
        @invariant !monom_isless(ht.monoms[lead_i], ht.monoms[lead_idx], ht.ord)
        !hashtable_monom_is_divisible(lead_i, lead_idx, ht) && continue

        # Add a new critical pair corresponding to Spoly(i, idx).
        pairset_resize_if_needed!(pairset, 1)
        ps[pairset.load + 1] = CriticalPair(Int32(idx), Int32(i), lead_i)
        degs[pairset.load + 1] = monom_totaldeg(ht.monoms[lead_i])
        pairset.load += 1

        basis.is_redundant[i] = true
    end
end

function basis_fill_data!(
    basis::Basis,
    ht::MonomialHashtable{M},
    exponents::Vector{Vector{M}},
    coeffs::Vector{Vector{T}}
) where {M, T}
    ngens = length(exponents)
    @inbounds for i in 1:ngens
        @invariant length(exponents[i]) == length(coeffs[i])
        hashtable_resize_if_needed!(ht, length(exponents[i]))

        nterms = length(coeffs[i])
        basis.coeffs[i] = coeffs[i]
        basis.monoms[i] = Vector{MonomId}(undef, nterms)
        poly = basis.monoms[i]
        @inbounds for j in 1:nterms
            poly[j] = hashtable_insert!(ht, exponents[i][j])
        end
    end

    basis.n_filled = ngens
end

function basis_discard_redundant_elements!(basis::Basis)
    j = 1
    @inbounds for i in 1:(basis.n_nonredundant)
        if !basis.is_redundant[basis.nonredundant_indices[i]]
            basis.divmasks[j] = basis.divmasks[i]
            basis.nonredundant_indices[j] = basis.nonredundant_indices[i]
            j += 1
        end
    end
    basis.n_nonredundant = j - 1
    @invariant basis.n_processed == basis.n_filled
    basis
end

function basis_standardize!(
    ring::PolyRing,
    basis::Basis,
    ht::MonomialHashtable,
    arithmetic::AbstractArithmetic,
    changematrix::Bool
)
    @inbounds for i in 1:(basis.n_nonredundant)
        idx = basis.nonredundant_indices[i]
        basis.nonredundant_indices[i] = i
        basis.is_redundant[i] = false
        basis.coeffs[i] = basis.coeffs[idx]
        basis.monoms[i] = basis.monoms[idx]
        if changematrix
            basis.changematrix[i] = basis.changematrix[idx]
        end
    end
    basis.n_processed = basis.n_filled = basis.n_nonredundant
    resize!(basis.coeffs, basis.n_processed)
    resize!(basis.monoms, basis.n_processed)
    resize!(basis.divmasks, basis.n_processed)
    resize!(basis.nonredundant_indices, basis.n_processed)
    resize!(basis.is_redundant, basis.n_processed)
    resize!(basis.changematrix, basis.n_processed)
    perm = sort_polys_by_lead_increasing!(basis, ht, changematrix, ord=ht.ord)
    basis_make_monic!(basis, arithmetic, changematrix)
    perm
end

function basis_get_monoms_by_identifiers(basis::Basis, ht::MonomialHashtable{M}) where {M <: Monom}
    monoms = Vector{Vector{M}}(undef, basis.n_nonredundant)
    @inbounds for i in 1:(basis.n_nonredundant)
        poly = basis.monoms[basis.nonredundant_indices[i]]
        monoms[i] = Vector{M}(undef, length(poly))
        for j in 1:length(poly)
            monoms[i][j] = ht.monoms[poly[j]]
        end
    end
    monoms
end

function basis_export_coeffs(basis::Basis{C}) where {C <: Coeff}
    coeffs = Vector{Vector{C}}(undef, basis.n_nonredundant)
    @inbounds for i in 1:(basis.n_nonredundant)
        coeffs[i] = basis.coeffs[basis.nonredundant_indices[i]]
    end
    coeffs
end

function basis_export_data(basis::Basis{C}, ht::MonomialHashtable{M}) where {M <: Monom, C <: Coeff}
    exps = basis_get_monoms_by_identifiers(basis, ht)
    coeffs = basis_export_coeffs(basis)
    exps, coeffs
end
