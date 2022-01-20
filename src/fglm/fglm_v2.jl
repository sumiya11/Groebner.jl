
#=
    The file contains the implementation of the FGLM algorithm,
    an algorithm for ordering conversion
=#

mutable struct NextMonomials
    # monomials to check
    monoms::Vector{Int}
    load::Int
end

Base.isempty(m::NextMonomials) = m.load == 0

function initialize_nextmonomials(ht)
    etmp = ht.exponents[1]
    etmp .= UInt16(0)
    vidx = insert_in_hash_table!(ht, etmp)
    monoms = Vector{Int}(undef, 2^64)
    monoms[1] = vidx
    load = 1
    NextMonomials(monoms, load)
end

function insertnexts!(
            m::NextMonomials, basis::Basis,
            ht::MonomialHashtable, monom::Int)
    etmp = ht.exponents[1]
    etmp .= UInt16(0)
    emonom = ht.exponents[monom]
    eprod = Vector{UInt16}(undef, ht.explen)
    for i in ht.nvars
        etmp[i] = UInt16(1)
        eprod = divisi
    end
end

function extractnext!(m::NextMonomials)
    # assuming m.monoms is sorted (reversed)
    monom = m.monoms[m.load]
    m.load -= 1
    monom
end

#------------------------------------------------------------------------------

function fglm_f4!(
            ring::PolyRing,
            basis::Basis,
            ht::MonomialHashtable)

    newbasis = initialize_basis(ring, basis.ntotal)
    nextmonoms = initialize_nextmonomials(ht)
    matrix = initialize_double_matrix()

    while !empty(nextmonoms)
        monom = extractnext!(nextmonoms)

        monom_nf, ht = normal_form_f4!(ring, basis, ht, monom)

        exists = linear_relation!(matrix, monom_nf)

        if exists
            add_generator!(newbasis, matrix)
        else
            insertnexts!(nextmonoms, newbasis, ht, monom)
        end
    end

    newbasis, ht
end

function fglm_f4(
            ring::PolyRing,
            basisexps::Vector{Vector{Vector{UInt16}}},
            basiscoeffs::Vector{Vector{UInt64}},
            rng::Rng,
            tablesize::Int=2^16)
    basis, ht = initialize_structures(ring, basisexps,
                                        basiscoeffs, rng, tablesize)

    basis, ht = fglm_f4!(ring, basis, ht, tobereduced)

    export_basis_data(basis, ht)
end
