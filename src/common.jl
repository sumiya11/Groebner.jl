
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

# Normal form of `h` with respect to `G`
function normal_form(h, G)
    i = 0
    while true
        if iszero(h)
            return h
        end
        i = 1
        while i <= length(G)
            mul = div(leading_monomial(h), leading_monomial(G[i]))
            if !iszero(mul)
                h -= leading_coefficient(h) * 1//leading_coefficient(G[i]) * mul * G[i]
                i = 1
                break
            end
            i = i + 1
        end
        if i > length(G)
            return h
        end
    end
end


function muls(f, g)
    lmi = leading_monomial(f)
    lmj = leading_monomial(g)
    lcm = lcm(lmi, lmj)
    mji = div(lcm, lmi)
    mij = div(lcm, lmj)
    return mji, mij
end

# generation of spolynomial of G[i] and G[j]
function spoly(f, g)
    mji, mij  = muls(f, g)
    h = 1//leading_coefficient(f) * mji * f - 1//leading_coefficient(g) * mij * g
    return h
end

function is_groebner(fs)
    for f in fs
        for g in fs
            if !iszero( normal_form(spoly(f, g)) )
                return false
            end
        end
    end
    return true
end



change_ordering(f, ordering) = change_ordering([f], ordering)

function change_ordering(fs::AbstractArray, ordering)
    R = parent(first(fs))
    Rord, _ = PolynomialRing(base_ring(R), string.(gens(R)), ordering=ordering)
    map(f -> change_base_ring(base_ring(R), f, parent=Rord), fs)
end
