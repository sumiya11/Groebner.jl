
mutable struct NextMonomials{Ord}
    # monomials to check
    monoms::Vector{MonomId}
    load::Int
    done::Dict{Int, Int}
    ord::Ord
end

Base.isempty(m::NextMonomials) = m.load == 0

function initialize_nextmonomials(ht::MonomialHashtable{M}, ord) where {M}
    zz = construct_const_monom(M, ht.nvars)
    vidx = insert_in_hashtable!(ht, zz)
    monoms = Vector{MonomId}(undef, 2^3)
    monoms[1] = vidx
    load = 1
    NextMonomials{typeof(ord)}(monoms, load, Dict{MonomId, Int}(vidx => 1), ord)
end

function insertnexts!(m::NextMonomials, ht::MonomialHashtable{M}, monom::MonomId) where {M}
    while m.load + ht.nvars >= length(m.monoms)
        resize!(m.monoms, length(m.monoms) * 2)
    end

    emonom = ht.monoms[monom]
    for i in 1:(ht.nvars)
        eprod = copy_monom(emonom)
        tmp = Vector{UInt}(undef, ht.nvars)
        edense = monom_to_dense_vector!(tmp, eprod)
        edense[i] += 1
        eprod = construct_monom(M, edense)
        vidx = insert_in_hashtable!(ht, eprod)

        if !haskey(m.done, vidx)
            m.load += 1
            m.monoms[m.load] = vidx
            m.done[vidx] = 1
        end
    end

    sort_monom_indices_decreasing!(m.monoms, m.load, ht, m.ord)
end

function nextmonomial!(m::NextMonomials)
    # assuming m.monoms is sorted (reversed)
    monom = m.monoms[m.load]
    m.load -= 1
    monom
end

function add_generator!(basis::Basis{C}, matrix, relation, ht, ord) where {C <: Coeff}
    rexps, rcoeffs, _ = extract_sparse_row(relation)

    for i in 1:length(rexps)
        rexps[i] = matrix.rightcolumn_to_monom[rexps[i]]
    end

    sort_term_indices_decreasing!(rexps, rcoeffs, ht, ord)

    resize_basis_if_needed!(basis, 1)
    basis.nprocessed += 1
    basis.nnonredundant += 1
    basis.nonredundant[basis.nnonredundant] = basis.nnonredundant
    basis.monoms[basis.nprocessed] = rexps
    basis.coeffs[basis.nprocessed] = rcoeffs
end

function divides_staircase(monom, staircase, ht)
    for m in staircase
        if is_monom_divisible(monom, m, ht)
            return true
        end
    end
    false
end

function fglm_f4!(
    ring::PolyRing,
    basis::Basis{C},
    ht::MonomialHashtable,
    ord::AbstractInternalOrdering,
    params
) where {C <: Coeff}
    newbasis = basis_initialize(ring, basis.nfilled, C)
    nextmonoms = initialize_nextmonomials(ht, ord)
    matrix = initialize_double_matrix(basis)
    staircase = MonomId[]

    while !isempty(nextmonoms)
        monom = nextmonomial!(nextmonoms)

        if divides_staircase(monom, staircase, ht)
            continue
        end

        tobereduced = basis_initialize(ring, [[monom]], [C[1]])
        tobereduced.nfilled = 1

        # compute normal form
        f4_normalform!(ring, basis, tobereduced, ht, params.arithmetic)

        # matrix left rows can express tobereduced?
        # reduces monom and tobereduced
        exists, relation =
            linear_relation!(ring, matrix, monom, tobereduced, ht, params.arithmetic)

        # if linear relation between basis elements exists
        if exists
            lead = ht.monoms[monom]
            add_generator!(newbasis, matrix, relation, ht, ord)
            push!(staircase, monom)
        else
            insertnexts!(nextmonoms, ht, monom)
        end
    end

    standardize_basis!(ring, newbasis, ht, ord, params.arithmetic)

    linbasis = extract_linear_basis(ring, matrix)
    newbasis, linbasis, ht
end

function fglm_f4(
    ring::PolyRing,
    basisexps::Vector{Vector{M}},
    basiscoeffs::Vector{Vector{C}},
    params
) where {M, C <: Coeff}
    basis, pairset, ht = initialize_structs(ring, basisexps, basiscoeffs, params)
    basis, linbasis, ht = fglm_f4!(ring, basis, ht, ring.ord, params)
    export_basis_data(basis, ht)
end

function extract_linear_basis(ring, matrix::DoubleMacaulayMatrix{C}) where {C}
    exps = Vector{Vector{MonomId}}(undef, matrix.nrrows)
    coeffs = Vector{Vector{C}}(undef, matrix.nrrows)

    for i in 1:(matrix.nrrows)
        exps[i] = matrix.rightrows[i]
        coeffs[i] = matrix.rightcoeffs[i]
        for j in 1:length(exps[i])
            exps[i][j] = matrix.rightcolumn_to_monom[exps[i][j]]
        end
    end

    linbasis = basis_initialize(ring, exps, coeffs)

    linbasis.nprocessed = length(exps)
    linbasis.nnonredundant = length(exps)
    linbasis.nonredundant = collect(1:(linbasis.nprocessed))

    linbasis
end
