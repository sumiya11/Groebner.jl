

import Nemo

function insert_nexts!(nextmonomials, monom)
    xs = gens(parent(monom))
    append!(nextmonomials, monom .* xs)
    sort!(unique!(nextmonomials))
end

function linear_relation(coeffs, vector)
    ground = parent(first(vector))
    m, n = size(coeffs)

    for k in 1:n

        for i in 1:m

            for j in 1:n

            end
        end

    end
end

function fglm(G)
    """


    """

    R = parent(first(G))
    xs = gens(R)
    ground = base_ring(R)

    newR, _ = PolynomialRing(ground, string.(xs), ordering=:lex)

    newG = MPoly[]
    staircase = MPoly[]
    MBasis = []

    nextmonomials = [ newR(1) ]
    i = 0

    while !isempty(nextmonomials)

        i += 1
        if i > 10
            break
        end

        m = pop!(nextmonomials)
        oldm = change_base_ring(ground, m, parent=R)

        if any(x -> divides(x, m)[1], staircase)
            continue
        end

        nf = normal_form(oldm, G)
        vector = collect(coefficients(nf))

        Mmatrix = [v for (x, v) in MBasis]

        println(Mmatrix)
        println(" / ")
        println(vector)
        println()

        λ, exists = linear_relation(Mmatrix, vector)
        println("$exists, $λ")

        if exists
            poly = oldm + sum(c * first(x) for (c, x) in zip(λ, MBasis))
            push!(newG, poly)
            push!(staircase, m)
        else
            push!(MBasis, (m, vector))
            insert_nexts!(nextmonomials, m)
        end

    end

    return newG
end
