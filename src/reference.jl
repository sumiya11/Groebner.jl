# Simple reference implementations for testing

#=
Given a Groebner basis G, constructs the reduced one.
Normalizes generators and sorts them by leading monomial.

Guarantees that If <f1..fn> equals <g1..gm> as ideals, then
    reducegb(groebner(f1..fn)) = reducegb(groebner(g1..gm))
=#
function _reducegb_reference(G)
    _reducegb_reference!(deepcopy(G))
end

function _reducegb_reference!(G)
    sort!(G, by=AbstractAlgebra.leading_monomial)
    for i in 1:length(G)
        G[i] = _normal_form_reference(G[i], G[1:end .!= i])
    end
    filter!(!iszero, G)
    scale = f -> AbstractAlgebra.map_coefficients(c -> c // AbstractAlgebra.leading_coefficient(f), f)
    # questionable over QQ
    map!(scale, G, G)
    sort!(G, by=AbstractAlgebra.leading_monomial)
end

#------------------------------------------------------------------------------

# Normal form of polynomial `h` with respect to ideal generators `G`
#
# The function is adapted from the eponymous function in
# https://www.mathematik.uni-kl.de/~ederc/teaching/2019/computeralgebra.html#news
function _normal_form_eder(h, G)
    i = 0
    while true
        if iszero(h)
            return h
        end
        i = 1
        while i <= length(G)
            mul = div(AbstractAlgebra.leading_monomial(h), AbstractAlgebra.leading_monomial(G[i]))
            if !iszero(mul)
                h -= AbstractAlgebra.leading_coefficient(h) * 1//AbstractAlgebra.leading_coefficient(G[i]) * mul * G[i]
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

# Normal form of polynomial `h` with respect to ideal generators `G`
function _normal_form_reference(h, G)
    R = parent(h)
    flag = false
    while true
        hprev = R(h)
        for t in AbstractAlgebra.terms(h)
            for g in G
                iszero(g) && continue
                is_divisible, _ = AbstractAlgebra.divides(t, AbstractAlgebra.leading_monomial(g))
                if is_divisible
                    mul = AbstractAlgebra.divexact(t, AbstractAlgebra.leading_monomial(g))
                    h -= mul * g
                    flag = true
                    break
                end
            end
            flag && break
        end
        !flag && break
        flag = false
        hprev == h && break
    end
    h
end

#------------------------------------------------------------------------------

function muls(f, g)
    lmi = AbstractAlgebra.leading_monomial(f)
    lmj = AbstractAlgebra.leading_monomial(g)
    lcm = AbstractAlgebra.lcm(lmi, lmj)
    mji = div(lcm, lmi)
    mij = div(lcm, lmj)
    mji, mij
end

# Generates an S-polynomial of f and g
#
# The function is adapted from the eponymous function in
# https://www.mathematik.uni-kl.de/~ederc/teaching/2019/computeralgebra.html#news
function _spoly_reference(f, g)
    mji, mij  = muls(f, g)
    h = 1//AbstractAlgebra.leading_coefficient(f) * mji * f - 1//AbstractAlgebra.leading_coefficient(g) * mij * g
    h
end

#------------------------------------------------------------------------------

# Checks if the given set of polynomials `fs` is a Groebner basis,
# i.e all spoly's reduce to zero.
# If `initial_gens` parameter is provided, also checks `initial_gens ⊆ fs` as ideals
function _isgroebner_reference(fs::Vector{AbstractAlgebra.Generic.MPoly{T}}; initial_gens=[]) where {T}
    sort!(fs, by=AbstractAlgebra.leading_monomial)
    for f in fs
        for g in fs
            if !iszero( _normal_form_reference(_spoly_reference(f, g), fs) )
                return false
            end
        end
    end
    if !isempty(initial_gens)
        return all(
            i -> _normal_form_reference( i, fs ) == 0,
            initial_gens
        )
    end
    true
end
