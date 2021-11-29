

import Nemo

function insert_nexts!(nextmonomials, monom)
    xs = sort(gens(parent(monom)))
    append!(nextmonomials, monom .* xs)
    sort!(unique!(nextmonomials))
end

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

function rowreduce!(Mvectors)
    m = length(Mvectors)
    n = length(first(Mvectors))
    pivots = Tuple{Int, Int}[]

    i, j = 1, 1
    while true
        ip, jp = findpivot(Mvectors, i, j)
        if ip == 0
            break
        end
        push!(pivots, (ip, jp))
        Mvectors[i], Mvectors[ip] = Mvectors[ip], Mvectors[i]
        Mvectors[i] = Mvectors[i] .* inv(Mvectors[i][jp])
        for k in 1:m
            if k == i
                continue
            end
            Mvectors[k] = Mvectors[k] .- Mvectors[k][jp] .* Mvectors[i]
        end
        i, j = ip+1, jp
    end
    Mvectors, pivots
end

function linear_relation!(Mvectors, vector)
    ground = parent(first(vector))
    Mreduced, pivots = rowreduce!(Mvectors)

    relation = zeros(ground, length(pivots))
    for (pivot, (ip, jp)) in enumerate(pivots)
        relation[pivot] = vector[jp]
        vector .-= vector[jp] .* Mreduced[ip]
    end

    return relation, iszero(vector)
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

function fglm(G)
    """


    """

    R = parent(first(G))
    xs = gens(R)
    ground = base_ring(R)

    newR, _ = PolynomialRing(ground, string.(xs), ordering=:lex)

    newG = []
    staircase = []
    MBasis = []

    monoms = []

    nextmonomials = [ newR(1) ]
    i = 0

    while !isempty(nextmonomials)

        i += 1
        if i > 10
            break
        end

        m = pop!(nextmonomials)
        oldm = change_base_ring(ground, m, parent=R)

        println("m = $m")
        if any(x -> divides(x, m)[1], staircase)
            continue
        end

        nf = R( GroebnerBases.normal_form(oldm, G) )
        union!(monoms, monomials(nf))

        println("nf = $nf, monoms = $monoms")
        vector = decompose_poly(nf, monoms)
        Mvectors = [decompose_poly(v, monoms) for (x, v) in MBasis]

        println(MBasis, " // ", Mvectors)
        println(vector)

        if !isempty(Mvectors)
            λ, exists = linear_relation!(Mvectors, vector)
        else
            exists = false
        end

        println("exists ?? $exists")

        if exists
            println("λ = $λ")

            poly = m + sum(c * first(x) for (c, x) in zip(λ, MBasis))
            push!(newG, poly)
            push!(staircase, m)
        else
            push!(MBasis, (m, nf))
            println(typeof(m), typeof(nf))
            union!(monoms, monomials(nf))
            insert_nexts!(nextmonomials, m)
        end

        println(newG, "\n", MBasis)
        println("--------------------")
    end

    return newG
end
