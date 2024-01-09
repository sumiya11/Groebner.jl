# This file is a part of Groebner.jl. License is GNU GPL v2.

# Enumerates monomials in the given monomial ordering
mutable struct MonomialEnumerator{Ord}
    # monomials to check
    monoms::Vector{MonomId}
    load::Int
    done::Dict{Int, Int}
    ord::Ord
end

Base.isempty(m::MonomialEnumerator) = m.load == 0

function enumerator_initialize(ht::MonomialHashtable{M}, ord) where {M}
    zz = monom_construct_const_monom(M, ht.nvars)
    vidx = hashtable_insert!(ht, zz)
    monoms = Vector{MonomId}(undef, 2^3)
    monoms[1] = vidx
    load = 1
    MonomialEnumerator{typeof(ord)}(monoms, load, Dict{MonomId, Int}(vidx => 1), ord)
end

function enumerator_produce_next_monomials!(
    m::MonomialEnumerator,
    ht::MonomialHashtable{M},
    monom::MonomId
) where {M}
    while m.load + ht.nvars >= length(m.monoms)
        resize!(m.monoms, length(m.monoms) * 2)
    end

    emonom = ht.monoms[monom]
    @inbounds for i in 1:(ht.nvars)
        eprod = monom_copy(emonom)
        tmp = Vector{UInt}(undef, ht.nvars)
        edense = monom_to_vector!(tmp, eprod)
        edense[i] += 1
        eprod = monom_construct_from_vector(M, edense)
        vidx = hashtable_insert!(ht, eprod)

        if !haskey(m.done, vidx)
            m.load += 1
            m.monoms[m.load] = vidx
            m.done[vidx] = 1
        end
    end

    sort_monom_indices_decreasing!(m.monoms, m.load, ht, m.ord)
end

function enumerator_next_monomial!(m::MonomialEnumerator)
    monom = m.monoms[m.load]
    m.load -= 1
    monom
end

function basis_add_generator!(basis::Basis, matrix, relation, ht, ord)
    rexps, rcoeffs, _ = extract_sparse_row(relation)

    for i in 1:length(rexps)
        rexps[i] = matrix.rightcolumn_to_monom[rexps[i]]
    end

    sort_term_indices_decreasing!(rexps, rcoeffs, ht, ord)

    basis_resize_if_needed!(basis, 1)
    basis.nprocessed += 1
    basis.nnonredundant += 1
    basis.nonredundant[basis.nnonredundant] = basis.nnonredundant
    basis.monoms[basis.nprocessed] = rexps
    basis.coeffs[basis.nprocessed] = rcoeffs
end

function staircase_divides_monom(monom, staircase, ht)
    for lead in staircase
        if hashtable_monom_is_divisible(monom, lead, ht)
            return true
        end
    end
    false
end

function fglm_main!(
    ring::PolyRing,
    basis::Basis{C},
    ht::MonomialHashtable{M},
    ord::AbstractInternalOrdering,
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @log level = -2 """
        Applying FGLM to convert a basis
        from ordering: $(ring.ord)
        to ordering:   $ord"""

    monom_enumerator = enumerator_initialize(ht, ord)
    new_basis = basis_initialize(ring, basis.nfilled, C)
    matrix = wide_matrix_initialize(basis)
    monom_staircase = Vector{MonomId}()

    # For "each" monomial, from the smallest to the largest
    while !isempty(monom_enumerator)
        # Get the next monomial
        monom = enumerator_next_monomial!(monom_enumerator)
        @log level = -4 "Processing monomial, id = $monom"

        # If the monomial is divisible by any of the elements of the existing
        # staircase, then throw the monomial away
        if staircase_divides_monom(monom, monom_staircase, ht)
            @log level = -4 "Monomial discarded, id = $monom"
            continue
        end

        # Compute the normal form of the monomial w.r.t. by constructing and
        # echelonizing the F4 matrix
        to_be_reduced = basis_initialize(ring, [[monom]], [C[1]])
        to_be_reduced.nfilled = 1
        f4_normalform!(ring, basis, to_be_reduced, ht, params.arithmetic)

        # Search for a linear relation between all the computed normal forms
        relation_exists, relation =
            find_linear_relation!(ring, matrix, monom, to_be_reduced, ht, params.arithmetic)

        if relation_exists
            # Add the monomial to the staircase, and add the pre-image of the
            # relation to the new basis
            push!(monom_staircase, monom)
            basis_add_generator!(new_basis, matrix, relation, ht, ord)
        end

        enumerator_produce_next_monomials!(monom_enumerator, ht, monom)
    end

    basis_standardize!(ring, new_basis, ht, ord, params.arithmetic)

    linbasis = extract_linear_basis(ring, matrix)
    new_basis, linbasis, ht
end
