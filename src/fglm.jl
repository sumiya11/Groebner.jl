

#=
    The file contains the implementation of the FGLM algorithm,
    an algorithm for ordering conversion
=#

import Nemo


# does not insert duplicates
function insertsorted!(array, item, rev=false)
    pos = searchsortedfirst(array, item, rev=rev)
    if !isempty(array) && pos <= length(array) && array[pos] == item
        return array
    end
    insert!(array, pos, item)
    return array
end

# TODO: revise this criterion
function insert_nexts!(nextmonomials, monom)
    xs = gens(parent(monom))
    for new_monom in monom .* xs
        insertsorted!(nextmonomials, new_monom)
    end
end

"""
    Represents the given polynomial as an element from Kⁿ
    vector space with respect to monomials occuring in monoms
"""
function decompose_poly(f, monoms)
    v = zeros(base_ring(f), length(monoms))
    for (j, c, m) in zip(1:length(f), coefficients(f), monomials(f))
        ind = findfirst(c -> c == m, monoms)
        v[ind] = c
    end
    v
end

"""
    Represents the given vector from Kⁿ as a polynomial
    in the given monomials
"""
function compose_poly(v, monoms)
    f = zero(parent(first(monoms)))
    for (i, m) in enumerate(monoms)
        f += m * v[i]
    end
    f
end

struct RrefBasis{C}
    rref_vectors::Vector{Vector{C}} # vectors stored in rref form
    monoms::Vector{MPoly{C}}   # monomials corresponding to columns of rref
    new_monoms::Vector{MPoly{C}} # monomials from the new basis
end

RrefBasis{C}() where {C} = RrefBasis(Vector{Vector{C}}(), Vector{MPoly{C}}(), Vector{MPoly{C}}())

function update_basis!(basis::RrefBasis, newpoly)
    ground = base_ring(newpoly)

    n = length(basis.monoms)
    for new_monom in monomials(newpoly)
        (new_monom in basis.monoms) && continue
        push!(basis.monoms, new_monom)
    end

    for vec in basis.rref_vectors
        append!(vec, zeros(ground, length(basis.monoms) - n))
    end
end

function solve_in_rref(rref_vectors, vector)
    coeffs = zeros(parent(first(vector)), length(rref_vectors))
    for (i, basis_vector) in enumerate(rref_vectors)
        pivot = findfirst(!iszero, basis_vector)
        iszero(vector[pivot]) && continue
        coeffs[i] = vector[pivot]
        vector .-= coeffs[i] .* basis_vector
    end
    exists = iszero(vector)
    coeffs, exists, vector
end

function linear_relation!(basis::RrefBasis, poly)
    vector = decompose_poly(poly, basis.monoms)
    coeffs, exists, vector = solve_in_rref(basis.rref_vectors, vector)
    coeffs, exists, vector
end

function update_rref!(basis::RrefBasis, new_monom, new_vector, coeffs)
    push!(basis.new_monoms, new_monom)
    pivot = findfirst(!iszero, new_vector)

    n_coeff = new_vector[pivot]
    new_vector .//= n_coeff


    for (i, c) in enumerate(coeffs)
        basis.new_monoms[end] -= c * basis.new_monoms[i]
    end
    basis.new_monoms[end] *= inv(n_coeff)


    es = basis.rref_vectors
    for i in 1:length(es)
        iszero(es[i][pivot]) && continue
        basis.new_monoms[i] -= es[i][pivot] * basis.new_monoms[end]
        es[i] .-= es[i][pivot] .* new_vector
    end
    push!(basis.rref_vectors, new_vector)

end

function linear_combination(basis::RrefBasis, coeffs)
    sum(c*x  for (c, x) in zip(coeffs, basis.new_monoms))
end

"""
    Converts the basis into the :lex ordering
"""
function fglm(G; linalg=:dense)

    origR  = parent(first(G))
    origxs = gens(origR)
    ground = base_ring(origR)
    T      = elem_type(origR)

    @assert ordering(origR) != :lex

    @info "converting from $(ordering(origR)) to :lex with fglm"

    newR, _ = PolynomialRing(ground, string.(origxs), ordering=:lex)

    newG      = T[]            #  lex
    staircase = T[]            #  orderless
    nextmonomials = T[newR(1)] # lex

    i = 0
    basis = RrefBasis{elem_type(ground)}()

    while !isempty(nextmonomials)
        i += 1
        if i > 10000
            @error "Something probably went wrong in FGLM"
            break
        end
        if i % 10 == 0
            @debug "fglm iter $i" size(newG)
        end

        m_lex = popfirst!(nextmonomials)
        m_deg = change_base_ring(ground, m_lex, parent=origR)

        # do we even need this ?
        if any(x -> is_term_divisible(m_lex, x), staircase)
            continue
        end

        # we convert to origR to ensure the output lives there
        nf = normal_form(m_deg, G)

        # pad vectors with zeros and add new monoms
        update_basis!(basis, nf)

        # try to find linear relation between vectors of basis
        # resulting to nf
        λ, exists, reduced = linear_relation!(basis, nf)

        if exists
            poly = m_lex - linear_combination(basis, λ)
            push!(newG, poly)
            push!(staircase, m_lex)
        else
            update_rref!(basis, m_lex, reduced, λ)
            insert_nexts!(nextmonomials, m_lex)
        end
    end

    return newG
end
