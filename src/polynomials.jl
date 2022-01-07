
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
struct 
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
            ht)

    e1 = ht[monom1]
    e2 = ht[monom2]

    N = length(e1)
    if ht.R.ord == :degrevlex
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
            ht)

    gcdexp = min.(ht[monom1], ht[monom2])

    if ht.R.ord == :degrevlex
        gcdexp[end] = 0
        for i in 1:length(gcdexp) - 1
            gcdexp[end] += gcdexp[i]
        end
    end

    h = add_monomial!(ht, gcdexp)

    h
end

# returns the lcm of given monomials in hashed representation
# and saves lcm to the given hash table
function monom_lcm(
            monom1::HashedMonomial,
            monom2::HashedMonomial,
            ht)

    lcmexp = max.(ht[monom1], ht[monom2])

    if ht.R.ord == :degrevlex
        lcmexp[end] = sum(lcmexp[1:end-1])
    end

    h = add_monomial!(ht, lcmexp)

    h
end

#------------------------------------------------------------------------------

# assesses if monom2 divides monom1
function is_monom_divisible(
            monom1::HashedMonomial,
            monom2::HashedMonomial,
            ht)

    e1 = ht[monom1]
    e2 = ht[monom2]
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
        ht) where {Tv}

    divisor = ht[monom1] .- ht[monom2]
    h = add_monomial!(ht, divisor)

    h
end

# divides term1 by term2 and returns the result in hash representation,
# saving the quotient in the given hashtable
function term_divide(
        term1::HashedTerm,
        term2::HashedTerm,
        ht) where {Tv}

    divisor = ht[term1] .- ht[term2]
    coeff   = term1.coeff // term2.coeff

    h = add_monomial!(ht, divisor)

    HashedTerm(h.e, coeff)
end

# divides monom1 by term2 and returns the result in hash representation,
# saving the quotient in the given hashtable
function term_divide(
        monom1::HashedMonomial,
        term2::HashedTerm,
        ht) where {Tv}

    divisor = ht[monom1] .- ht[HashedMonomial(term2.e)]
    coeff   = 1 // term2.coeff

    h = add_monomial!(ht, divisor)

    HashedTerm(h.e, coeff)
end

#------------------------------------------------------------------------------

# computes the product monom*f in hash representation
# and stores intermediate monomials in hash table
function monom_poly_product(
        monom::HashedMonomial,
        f::HashedPolynomial,
        ht)

    hm = ht[monom]
    exps = [
        add_monomial!(ht, hm .+ ht[e])
        for e in f.monoms
    ]

    #=
    exps = [
        add_product!(ht, monom, e)
        for e in f.monoms
    ]=#

    HashedPolynomial(exps, f.coeffs, f.parent)
end

# computes the product term*f in hash representation
# and stores intermediate monomials in hash table
function term_poly_product(
        term::HashedTerm,
        f::HashedPolynomial,
        ht)

    coeffs = term.coeff .* f.coeffs

    #=
    for e in f.monoms
        ee = ht[HashedMonomial(term.e)] .+ ht[e]
        if haskey(ht.expmap, hash(ee))
            @warn "hit prod!"
        end
    end
    =#

    exps = [
        add_monomial!(ht, ht[HashedMonomial(term.e)] .+ ht[e])
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
        ht)

    A1, A2 = ht[monom1], ht[monom2]
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
                ht) where {Tv, Ti}

    len  = nnz(vector)
    raw_explen = length(ht[first(monom_basis)])

    cfs  = Vector{Tv}(undef, len)
    exps = Vector{HashedMonomial}(undef, len)

    for (i, midx, val) in zip(1:len, nonzeroinds(vector), nonzeros(vector))
        cfs[i]  = val
        exps[i] = monom_basis[midx]
    end

    HashedPolynomial(exps, cfs, R)
end

#------------------------------------------------------------------------------

function convert_to_hash_repr(orig_poly::MPoly{Tv}, ht) where {Tv}
    n = length(orig_poly)
    imonoms = Vector{HashedMonomial}(undef, n)
    icoeffs = Vector{Tv}(undef, n)
    for (i, c) in zip(1:n, coefficients(orig_poly))
        h = add_monomial!(ht, orig_poly.exps[:, i])
        imonoms[i] = h
        icoeffs[i] = c
    end
    HashedPolynomial(imonoms, icoeffs, parent(orig_poly))
end

function convert_to_hash_repr(polys::Vector{MPoly{Tv}}, ht) where {Tv}
    ans = Vector{}(undef, length(polys))
    for i in 1:length(polys)
        ans[i] = convert_to_hash_repr(polys[i], ht)
    end
    ans
end

#------------------------------------------------------------------------------

function convert_to_original_repr(orig_poly::, ht) where {Tv}
    n = length(orig_poly)
    iexps = MacaulayMatrix{UInt}(undef, length(ht[first(orig_poly.monoms)]), n)
    icoeffs = Vector{Tv}(undef, n)
    for (i, c, m) in zip(1:n, coefficients(orig_poly), monomials(orig_poly))
        iexps[:, i] .= ht[m]
        icoeffs[i] = c
    end
    MPoly{Tv}(parent(orig_poly), icoeffs, iexps)
end

function convert_to_original_repr(polys::Vector{}, ht) where {Tv}
    ans = Vector{MPoly{Tv}}(undef, length(polys))
    for i in 1:length(polys)
        ans[i] = convert_to_original_repr(polys[i], ht)
    end
    ans
end

#------------------------------------------------------------------------------
