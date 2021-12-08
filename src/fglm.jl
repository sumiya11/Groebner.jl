

#=
    The file contains the implementation of the FGLM algorithm,
    an algorithm for ordering conversion
=#

import Nemo


function findpivot(Mvectors, row, col)
    m = length(Mvectors)
    n = length(first(Mvectors))
    i, j = row, col
    while j <= n
        while i <= m
            if Mvectors[i][j] != 0
                return (i, j)
            end
            i += 1
        end
        i = row
        j += 1
    end
    return (0, 0)
end

function rowreduce!(Mvectors, vector)
    m, n = length(Mvectors), length(first(Mvectors))
    pivots = Tuple{Int, Int}[]

    i, j = 1, 1
    while true
        ip, jp = findpivot(Mvectors, i, j)
        if ip == 0
            break
        end

        Mvectors[i], Mvectors[ip] = Mvectors[ip], Mvectors[i]
        vector[i], vector[ip] = vector[ip], vector[i]

        vector[i] = vector[i] .// Mvectors[i][jp]
        Mvectors[i] = Mvectors[i] .// Mvectors[i][jp]

        push!(pivots, (i, jp))

        for k in i+1:m
            if k == i
                continue
            end
            vector[k] = vector[k] - vector[i] * Mvectors[k][jp]
            Mvectors[k] = Mvectors[k] .- Mvectors[i] .* Mvectors[k][jp]
        end
        i, j = i + 1, jp
    end

    for (i, j) in pivots
        for k in 1:i-1
            vector[k] = vector[k] - vector[i] * Mvectors[k][j]
            Mvectors[k] = Mvectors[k] .- Mvectors[i] .* Mvectors[k][j]
        end
    end

    Mvectors, vector, pivots
end

function linear_relation!(Mvectors, vector)
    ground = parent(first(vector))

    Mmatrix = [hcat(Mvectors...)[i, :] for i in 1:length(vector)]
    Mreduced, vector, pivots = rowreduce!(Mmatrix, vector)

    λ = zeros(ground, length(Mvectors))

    for (ip, jp) in pivots
        λ[jp] = vector[ip]
    end

    return λ, all(iszero, vector[length(pivots)+1:end])
end


# TODO: revise this criterion
function insert_nexts!(nextmonomials, monom)
    xs = gens(parent(monom))
    append!(nextmonomials, monom .* xs)
    sort!(unique!(nextmonomials), rev=false)
end

function insert_nexts_2!(nextmonomials, monom, state=(1, 0))
    xs = sort(gens(parent(monom)))[state[1]]
    push!(nextmonomials, xs^state[2])
    sort!(unique!(nextmonomials), rev=false)
end


"""
    Represents the given polynomial as an element from Kⁿ
    vector space with respect to monomials occuring in monoms
"""
function decompose_poly(f, monoms)
    v = zeros(base_ring(f), length(monoms))
    for (i, m) in enumerate(monoms)
        v[i] = coeff(f, m)
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


"""
    Converts the basis into the :lex ordering
"""
function fglm(G)
    origR = parent(first(G))
    origxs = gens(origR)
    ground = base_ring(origR)

    @assert ordering(origR) != :lex

    @info "converting from $(ordering(origR)) to :lex with fglm"

    newR, _ = PolynomialRing(ground, string.(origxs), ordering=:lex)

    newG = []      #  lex
    staircase = [] #  orderless
    MBasis = []    #  <a, b>   a = nf wrt newBasis, b = nf(a) wrt oldBasis
    nextmonomials = [ newR(1) ] # lex

    i = 0
    monoms = []    # orderless (?)

    while !isempty(nextmonomials)
        i += 1
        if i > 1000
            @error "Something is probably wrong in FGLM"
            break
        end

        m_lex = popfirst!(nextmonomials)
        m_deg = change_base_ring(ground, m_lex, parent=origR)

        if any(x -> divides(m_lex, x)[1], staircase)
            continue
        end

        # we convert to origR to ensure the output lives there
        nf = origR( normal_form(m_deg, G) )
        union!(monoms, monomials(nf))
        sort!(monoms, rev=false)

        vector = decompose_poly(nf, monoms)
        Mvectors = [decompose_poly(v, monoms) for (x, v) in MBasis]

        if !isempty(Mvectors)
            λ, exists = linear_relation!(Mvectors, vector)
        else
            exists = false
        end

        if exists
            poly = m_lex - sum(c * first(x) for (c, x) in zip(λ, MBasis))
            push!(newG, poly)
            push!(staircase, m_lex)

        else
            push!(MBasis, (m_lex, nf))
            # # println(typeof(m_lex), typeof(nf))
            union!(monoms, monomials(nf))

            insert_nexts!(nextmonomials, m_lex)
        end
    end

    return newG
end
