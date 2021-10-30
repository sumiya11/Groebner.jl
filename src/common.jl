
# The file contains definitions of some common functions


# Constructs reduced groebner basis given that G is itself a groebner basis,
# strightforward approach :D
# Note: groebner bases generating same ideal meet same reduced form
function reducegb(G)
    GG = deepcopy(G)
    n = length(G)
    updated = true
    while updated
        updated = false
        for (i, f) in enumerate(GG)
            for m in collect(terms(f))
                for (j, g) in enumerate(GG[1:end .!= i])
                    flag, t = divides(m, leading_term(g))
                    if flag
                        GG[i] -= t*g
                        updated = true
                        break
                    end
                end
            end
        end
    end
    GG = filter!(!iszero, GG)
    GG = map(f -> map_coefficients(c -> c // leading_coefficient(f), f), GG)
    sort(GG, by=leading_monomial)
end
