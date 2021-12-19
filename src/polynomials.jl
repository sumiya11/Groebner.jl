
import Base: ==

#------------------------------------------------------------------------------

# A hash-based representation for monomial. Provides usual monomial operations
# interface. The corresponding exponent should be stored in the current
# ::MonomialHashtable object
struct HashedMonomial
    e::UInt # hashed exponent
end

# we will compare them fast !
==(m1::HashedMonomial, m2::HashedMonomial) = (m1.e == m2.e)
Base.hash(m1::HashedMonomial) = m1.e

#------------------------------------------------------------------------------

# A hash-based representation for a term
struct HashedTerm{Tv}
    e::UInt   #  hashed exponent
    coeff::Tv # coefficient
end

# A hashed-based representation for a polynomial.
# Internally stores hashed monomials, coefficients, and base ring
struct HashedPolynomial{Tv}
    monoms::Vector{HashedMonomial}  # hashed monoms, sorted wrt ordering(parent)
    coeffs::Vector{Tv}              # coeffs, -//-
    parent::MPolyRing{Tv}
end

#------------------------------------------------------------------------------

function AbstractAlgebra.leading_monomial(f::HashedPolynomial)
    return first(f.monoms)
end

function AbstractAlgebra.leading_term(f::HashedPolynomial)
    return HashedTerm(first(f.monoms).e, first(f.coeffs))
end

# semantic of function is changed
function AbstractAlgebra.monomials(f::HashedPolynomial)
    return f.monoms
end

# semantic of function is changed
function AbstractAlgebra.coefficients(f::HashedPolynomial)
    return f.coeffs
end

function AbstractAlgebra.base_ring(f::HashedPolynomial)
    return base_ring(f.parent)
end

function Base.length(f::HashedPolynomial)
    return length(f.monoms)
end

function Base.parent(f::HashedPolynomial)
    return f.parent
end

#------------------------------------------------------------------------------

# checks if given monomials are coprime
function is_monom_gcd_constant(
            monom1::HashedMonomial,
            monom2::HashedMonomial,
            HT)

    e1 = HT[monom1]
    e2 = HT[monom2]

    N = length(e1)
    if HT.R.ord == :degrevlex
        N -= 1
    end

    for i in 1:N
        if min(e1[i], e2[i]) != 0
            return false
        end
    end

    return true
end

# returns the gcd of given monomials in hashed representation
# and saves gcd to the given hash table
function monom_gcd(
            monom1::HashedMonomial,
            monom2::HashedMonomial,
            HT)

    gcdexp = min.(HT[monom1], HT[monom2])

    if HT.R.ord == :degrevlex
        gcdexp[end] = 0
        for i in 1:length(gcdexp) - 1
            gcdexp[end] += gcdexp[i]
        end
    end

    h = add_monomial!(HT, gcdexp)

    h
end

# returns the lcm of given monomials in hashed representation
# and saves lcm to the given hash table
function monom_lcm(
            monom1::HashedMonomial,
            monom2::HashedMonomial,
            HT)

    lcmexp = max.(HT[monom1], HT[monom2])

    if HT.R.ord == :degrevlex
        lcmexp[end] = sum(lcmexp[1:end-1])
    end

    h = add_monomial!(HT, lcmexp)

    h
end

#------------------------------------------------------------------------------

# assesses if monom2 divides monom1
function is_monom_divisible(
            monom1::HashedMonomial,
            monom2::HashedMonomial,
            HT)

    e1 = HT[monom1]
    e2 = HT[monom2]
    for i in 1:length(e1)
        if e1[i] < e2[i]
            return false
        end
    end
    return true
end

# divides monom1 by monom2 and returns the result in hash representation,
# saving the quotient in the given hashtable
function monom_divide(
        monom1::HashedMonomial,
        monom2::HashedMonomial,
        HT) where {Tv}

    divisor = HT[monom1] .- HT[monom2]
    h = add_monomial!(HT, divisor)

    h
end

# divides term1 by term2 and returns the result in hash representation,
# saving the quotient in the given hashtable
function term_divide(
        term1::HashedTerm,
        term2::HashedTerm,
        HT) where {Tv}

    divisor = HT[term1] .- HT[term2]
    coeff   = term1.coeff // term2.coeff

    h = add_monomial!(HT, divisor)

    HashedTerm(h.e, coeff)
end

# divides monom1 by term2 and returns the result in hash representation,
# saving the quotient in the given hashtable
function term_divide(
        monom1::HashedMonomial,
        term2::HashedTerm,
        HT) where {Tv}

    divisor = HT[monom1] .- HT[HashedMonomial(term2.e)]
    coeff   = 1 // term2.coeff

    h = add_monomial!(HT, divisor)

    HashedTerm(h.e, coeff)
end

#------------------------------------------------------------------------------

# computes the product monom*f in hash representation
# and stores intermediate monomials in hash table
function monom_poly_product(
        monom::HashedMonomial,
        f::HashedPolynomial,
        HT)

    hm = HT[monom]
    exps = [
        add_monomial!(HT, hm .+ HT[e])
        for e in f.monoms
    ]

    #=
    exps = [
        add_product!(HT, monom, e)
        for e in f.monoms
    ]=#

    HashedPolynomial(exps, f.coeffs, f.parent)
end

# computes the product term*f in hash representation
# and stores intermediate monomials in hash table
function term_poly_product(
        term::HashedTerm,
        f::HashedPolynomial,
        HT)

    coeffs = term.coeff .* f.coeffs

    #=
    for e in f.monoms
        ee = HT[HashedMonomial(term.e)] .+ HT[e]
        if haskey(HT.expmap, hash(ee))
            @warn "hit prod!"
        end
    end
    =#

    exps = [
        add_monomial!(HT, HT[HashedMonomial(term.e)] .+ HT[e])
        for e in f.monoms
    ]

    HashedPolynomial(exps, coeffs, f.parent)
end

#------------------------------------------------------------------------------

# True if monom1 < monom2 in ordering imposed by ordering(R)
# False otherwise
function monom_is_less(
        monom1::HashedMonomial,
        monom2::HashedMonomial,
        R,
        HT)

    A1, A2 = HT[monom1], HT[monom2]
    if R.ord == :degrevlex
      N = nvars(R) + 1
      if A1[N] < A2[N]
         return true
     elseif A1[N] > A2[N]
         return false
      end
      for k = N-1:-1:1
         if A1[k] > A2[k]
            return true
        elseif A1[k] < A2[k]
            return false
         end
      end
   else
      N = nvars(R)
      for k = N:-1:1
         if A1[k] < A2[k]
            return true
        elseif A1[k] > A2[k]
            return false
         end
      end
   end
   return false
end

#------------------------------------------------------------------------------

# converts sparse polynomial representation wrt some monomial basis monom_basis
# to hash format. Therefore all exponents are guaranteed to be hashed before
function sparsevector_to_hashpoly(
                vector::SparseVectorAA{Tv, Ti},
                monom_basis,
                R,
                HT) where {Tv, Ti}

    len  = nnz(vector)
    raw_explen = length(HT[first(monom_basis)])

    cfs  = Vector{Tv}(undef, len)
    exps = Vector{HashedMonomial}(undef, len)

    for (i, midx, val) in zip(1:len, nonzeroinds(vector), nonzeros(vector))
        cfs[i]  = val
        exps[i] = monom_basis[midx]
    end

    HashedPolynomial(exps, cfs, R)
end

#------------------------------------------------------------------------------

function convert_to_hash_repr(orig_poly::MPoly{Tv}, HT) where {Tv}
    n = length(orig_poly)
    imonoms = Vector{HashedMonomial}(undef, n)
    icoeffs = Vector{Tv}(undef, n)
    for (i, c) in zip(1:n, coefficients(orig_poly))
        h = add_monomial!(HT, orig_poly.exps[:, i])
        imonoms[i] = h
        icoeffs[i] = c
    end
    HashedPolynomial(imonoms, icoeffs, parent(orig_poly))
end

function convert_to_hash_repr(polys::Vector{MPoly{Tv}}, HT) where {Tv}
    ans = Vector{HashedPolynomial{Tv}}(undef, length(polys))
    for i in 1:length(polys)
        ans[i] = convert_to_hash_repr(polys[i], HT)
    end
    ans
end

#------------------------------------------------------------------------------

function convert_to_original_repr(orig_poly::HashedPolynomial{Tv}, HT) where {Tv}
    n = length(orig_poly)
    iexps = Matrix{UInt}(undef, length(HT[first(orig_poly.monoms)]), n)
    icoeffs = Vector{Tv}(undef, n)
    for (i, c, m) in zip(1:n, coefficients(orig_poly), monomials(orig_poly))
        iexps[:, i] .= HT[m]
        icoeffs[i] = c
    end
    MPoly{Tv}(parent(orig_poly), icoeffs, iexps)
end

function convert_to_original_repr(polys::Vector{HashedPolynomial{Tv}}, HT) where {Tv}
    ans = Vector{MPoly{Tv}}(undef, length(polys))
    for i in 1:length(polys)
        ans[i] = convert_to_original_repr(polys[i], HT)
    end
    ans
end

#------------------------------------------------------------------------------
