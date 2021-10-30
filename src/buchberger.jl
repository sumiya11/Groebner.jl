"""
    The file is taken from the Coumputer Algebra course by Christian Eder,
    https://www.mathematik.uni-kl.de/~ederc/teaching/2019/computeralgebra.html#news
    and adapted to avoid C dependencies by Alexander Demin
"""

import AbstractAlgebra

zero_reductions = 0
non_zero_reductions = 0
prod_criteria = 0

# structure for pairset
mutable struct Pairset
    len::Int64
    pairs::Array{Tuple{Int64,Int64}, 1}
end

# criteria for useless pairs
function product_criterion(pair::Tuple{Int64,Int64}, G)
    lm1 = AbstractAlgebra.leading_monomial(G[pair[1]])
    lm2 = AbstractAlgebra.leading_monomial(G[pair[2]])
    if AbstractAlgebra.lcm(lm1, lm2) == lm1 * lm2
        return true
    end
    return false
end

# initialize pairset
function init_pairset(idx::Int64)
    len = 0
    for k in 1:idx-1
        len +=  k
    end
    p   = Array{Tuple{Int64,Int64},1}(undef, len)
    P = Pairset(len, p)
    k = 1
    for i = 2:idx
        for j = 1:i-1
            P.pairs[k]  = (i,j)
            k = k+1
        end
    end
    return P
end

# add new pairs to pairset
function new_pairs(P::Pairset, idx::Int64)
    k     = P.len + 1
    P.len = P.len + idx - 1
    resize!(P.pairs, P.len)
    for i in 1:idx-1
        P.pairs[k]  = (idx,i)
        k = k+1
    end
end

# multiplier for reduction
function multiplier_exp(h, g)
    x   = AbstractAlgebra.gens(parent(h))
    eh  = AbstractAlgebra.lead_exponent(h)
    eg  = AbstractAlgebra.lead_exponent(g)

    mul = one(parent(h))

    for k in 1:length(eh)
        diff = eh[k] - eg[k]
        if diff < 0
            return 0
        end
        if diff > 0
            mul *=  x[k]^diff
            continue
        end
    end
    return mul
end

# multipliers for spolynomials
function muls_exp(i::Int64, j::Int64, G)
    x   = AbstractAlgebra.gens(parent(G[1]))
    ei  = AbstractAlgebra.lead_exponent(G[i])
    ej  = AbstractAlgebra.lead_exponent(G[j])

    mji = one(parent(G[1]))
    mij = one(parent(G[1]))
    for k in 1:length(ei)
        eij = ej[k] - ei[k]
        if eij > 0
            mji *=  x[k]^eij
            continue
        end
        if eij < 0
            mij *= x[k]^(-eij)
        end
    end
    return mji, mij
end

function muls(i::Int64, j::Int64, G)
    lmi = AbstractAlgebra.leading_monomial(G[i])
    lmj = AbstractAlgebra.leading_monomial(G[j])
    lcm = AbstractAlgebra.lcm(lmi, lmj)
    mji = AbstractAlgebra.div(lcm, lmi)
    mij = AbstractAlgebra.div(lcm, lmj)
    return mji, mij
end

# generation of spolynomial of G[i] and G[j]
function spoly(pair::Tuple{Int64,Int64}, G)
    i   = pair[1]
    j   = pair[2]

    mji, mij  = muls(i, j, G)
    h = 1//AbstractAlgebra.leading_coefficient(G[i]) * mji * G[i] - 1//AbstractAlgebra.leading_coefficient(G[j]) * mij * G[j]

    return h
end

# compute normal form resp. standard expression of h w.r.t. G
function normal_form(h, G)
    i = 0
    while true
        if h == 0
            global zero_reductions += 1
            return 0
        end
        i = 1
        while i <= length(G)
            mul = AbstractAlgebra.div(AbstractAlgebra.leading_monomial(h), AbstractAlgebra.leading_monomial(G[i]))
            # mul = multiplier_exp(h, G[i])
            if mul != 0
                h = h - AbstractAlgebra.leading_coefficient(h) * 1//AbstractAlgebra.leading_coefficient(G[i]) * mul * G[i]
                i = 1
                break
            end
            i = i + 1
        end
        if i > length(G)
            global non_zero_reductions += 1
            return h
        end
    end
end

function bba(id::Array{T}) where {T}
    global zero_reductions = 0
    global non_zero_reductions = 0
    global prod_criteria  = 0
    # get number of generators in ideal
    R = AbstractAlgebra.parent(first(id))
    l = AbstractAlgebra.nvars(R)
    # array for storing gb
    G = Array{T,1}(undef, l)
    # put all initial generators in G
    for i in 1:l
        G[i] = id[i]
    end

    # generate list of tuples of indices for s-polynomials
    P = init_pairset(l)

    # as long as the pairset is not empty
    while P.len != 0
        # always take the last element and adjust length of P
        pair  = P.pairs[P.len]
        P.len = P.len - 1
        if product_criterion(pair, G)
            global prod_criteria += 1
            continue
        end
        h0     = spoly(pair, G)
        h     = normal_form(h0, G)

        # @info "step" h0 h

        # if Buchberger's criterion is not fulfilled
        if h != 0
            new_pairs(P, length(G)+1)
            push!(G, h)
        end
    end

    # return gb with setting isGB flag set

    println("#zero reductions     ", zero_reductions)
    println("#non-zero reductions ", non_zero_reductions)
    println("#product criteria    ", prod_criteria)
    return sort(G, by=leading_monomial)
end


#=
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
fs = [ x + y + z, x*y + x*z + y*z, x*y*z - 1 ]

gb = bba(fs)

println(gb)
=#
