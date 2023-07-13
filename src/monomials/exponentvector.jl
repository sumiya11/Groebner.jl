# ExponentVector{T} is a dense vector of integers of type T that implements the
# monomial interface.

# ExponentVector stores the total degree inline with the partial degrees.
# For example, x^2 y z^5 is stored as a dynamic vector [8, 2, 1, 5].

# ExponentVector is just an alias for a dynamic vector of integers.
const ExponentVector{T} = Vector{T} where {T <: Integer}

# Checks whether ExponentVector{T} provides an efficient comparator function
# implementation for the given monomial ordering of type `O`
function is_supported_ordering(::Type{ExponentVector{T}}, ::O) where {T, O}
    # ExponentVector{T} supports efficient implementation of all monomial
    # orderings
    true
end

# Checks if there is a risk of exponent overflow. If overflow if possible,
# throws a MonomialDegreeOverflow.
function _monom_overflow_check(e::ExponentVector{T}) where {T}
    # ExponentVector overflows if its total degree overflows
    _monom_overflow_check(totaldeg(e), T)
end

# The maximum number of variables that this monomial implementation can
# potentially store. 
max_vars_in_monom(::Type{ExponentVector{T}}) where {T} = 2^32
max_vars_in_monom(p::ExponentVector{T}) where {T} = max_vars_in_monom(typeof(p))

# The total degree of a monomial
totaldeg(pv::ExponentVector) = @inbounds pv[1]

# The type of an entry of a ExponentVector{T}. Note that this is not necessarily
# equal to T.
# NOTE: this may be a bit awkward
entrytype(::Type{ExponentVector{T}}) where {T} = MonomHash
entrytype(p::ExponentVector{T}) where {T} = entrytype(typeof(p))

# Constructs a constant monomial of type ExponentVector{T} with the max_vars_in_monom for
# n variables.
construct_const_monom(::Type{ExponentVector{T}}, n::Integer) where {T} = zeros(T, n + 1)

# Constructs a monomial of type ExponentVector{T} with the partial degrees taken
# from the vector `ev`
function construct_monom(::Type{ExponentVector{T}}, ev::Vector{U}) where {T, U}
    v = Vector{T}(undef, length(ev) + 1)
    @inbounds v[1] = sum(ev)
    @inbounds v[2:end] .= ev
    ExponentVector{T}(v)
end

# Returns a dense vector of partial degrees that correspond to the monomial `pv`.
function monom_to_dense_vector!(tmp::Vector{M}, pv::ExponentVector{T}) where {M, T}
    @assert length(tmp) == length(pv) - 1
    @inbounds tmp[1:end] = pv[2:end]
    tmp
end

# Returns a monomial that can be used to compute hashes of monomials of type
# ExponentVector{T}
function construct_hash_vector(::Type{ExponentVector{T}}, n::Integer) where {T}
    rand(MonomHash, n + 1)
end

# Computes the hash of `x` with the given hash vector `b`
function monom_hash(x::ExponentVector{T}, b::ExponentVector{MH}) where {T, MH}
    h = zero(MH)
    @inbounds for i in eachindex(x, b)
        h += MH(x[i]) * b[i]
    end
    mod(h, MonomHash)
end

###
# Monomial comparator functions. See monoms/orderings.jl for details.

# DegRevLex monomial comparison
function monom_isless(ea::ExponentVector, eb::ExponentVector, ::DegRevLex)
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    if @inbounds ea[1] < eb[1]
        return true
    elseif @inbounds ea[1] != eb[1]
        return false
    end
    i = length(ea)
    @inbounds while i > 2 && ea[i] == eb[i]
        i -= 1
    end
    @inbounds if ea[i] <= eb[i]
        return false
    else
        return true
    end
end

# DegRevLex monomial comparison, but on a given range of variables [lo:hi]. The
# flag `hasdegree` indicates whether the total degree is stored in ea[lo] and
# eb[lo], respectively.
function monom_isless(
    ea::ExponentVector,
    eb::ExponentVector,
    ::DegRevLex,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    @invariant length(ea) == length(eb)
    @invariant hi >= lo
    @inbounds da = hasdegree ? ea[lo] : sum(view(ea, lo:hi))
    @inbounds db = hasdegree ? eb[lo] : sum(view(eb, lo:hi))
    if da < db
        return true
    elseif da != db
        return false
    end
    # after the while loop terminates, 
    # i is guaranteed to be in [lo + 1, hi]
    i = hi
    @inbounds while i > lo && ea[i] == eb[i]
        i -= 1
    end
    @inbounds if ea[lo] <= eb[lo]
        return false
    else
        return true
    end
end

# DegLex exponent vector comparison
function monom_isless(ea::ExponentVector, eb::ExponentVector, ::DegLex)
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    if @inbounds ea[1] < eb[1]
        return true
    elseif @inbounds ea[1] != eb[1]
        return false
    end
    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end

# DegLex monomial comparison, but on a given range of variables [lo:hi]
function monom_isless(
    ea::ExponentVector,
    eb::ExponentVector,
    ::DegLex,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    @invariant length(ea) == length(eb)
    @invariant hi >= lo
    @inbounds da = hasdegree ? ea[lo] : sum(view(ea, lo:hi))
    @inbounds db = hasdegree ? eb[lo] : sum(view(eb, lo:hi))
    if da < db
        return true
    elseif da != db
        return false
    end
    i = lo + hasdegree
    @inbounds while i < hi && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end

# Lex monomial comparison
function monom_isless(ea::ExponentVector, eb::ExponentVector, ::Lex)
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end

# Lex monomial comparison, but on a given range of variables [lo:hi].
function monom_isless(
    ea::ExponentVector,
    eb::ExponentVector,
    ::Lex,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    @invariant length(ea) == length(eb)
    @invariant hi >= lo
    i = lo + hasdegree
    @inbounds while i < hi && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end

#
# Weighted monomial comparison
function monom_isless(ea::ExponentVector, eb::ExponentVector, w::WeightedOrdering)
    monom_isless(ea, eb, w, 1, length(ea), true)
end
# Weighted monomial comparison, but on a given range of variables [lo:hi].
function monom_isless(
    ea::ExponentVector,
    eb::ExponentVector,
    w::WeightedOrdering,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    @invariant length(ea) == length(eb)
    @invariant hi >= lo
    weights = w.weights
    common_type = promote_type(eltype(weights), eltype(ea))
    sa, sb = zero(common_type), zero(common_type)
    j = 1
    @inbounds for i in (lo + hasdegree):hi
        # TODO: check for overflow!
        sa += weights[j] * ea[i]
        sb += weights[j] * eb[i]
        j += 1
    end
    if sa < sb
        true
    elseif sa == sb
        false
    else
        false
    end
end

# Block exponent vector comparison.
function monom_isless(ea::ExponentVector, eb::ExponentVector, b::BlockOrdering)
    monom_isless(ea, eb, b, 1, length(ea), true)
end
# Block, but on a given range of variables [lo:hi].
function monom_isless(
    ea::ExponentVector,
    eb::ExponentVector,
    b::BlockOrdering,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    # here, it is assumed that lo:hi equals the union of r1 and r2
    r1, r2 = b.r1, b.r2
    if monom_isless(ea, eb, b.ord1, first(r1) + hasdegree, last(r1) + hasdegree, false)
        true
    elseif monom_isless(eb, ea, b.ord1, first(r1) + hasdegree, last(r1) + hasdegree, false)
        false
    else
        monom_isless(ea, eb, b.ord2, first(r2) + hasdegree, last(r2) + hasdegree, false)
    end
end

# Matrix exponent vector comparison.
function monom_isless(ea::ExponentVector, eb::ExponentVector, m::MatrixOrdering)
    monom_isless(ea, eb, m, 1, length(ea), true)
end
# Matrix, but on a given range of variables [lo:hi].
function monom_isless(
    ea::ExponentVector,
    eb::ExponentVector,
    m::MatrixOrdering,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    rows = m.rows
    @inbounds common_type = promote_type(eltype(rows[1]), eltype(ea))
    @inbounds for i in 1:length(rows)
        sa, sb = zero(common_type), zero(common_type)
        k = 1
        for j in (lo + hasdegree):hi
            sa += rows[i][k] * common_type(ea[j])
            sb += rows[i][k] * common_type(eb[j])
            k += 1
        end
        if sa < sb
            return true
        elseif sa > sb
            return false
        end
    end
    false
end

###
# Monomial-Monomial arithmetic

# Returns the lcm of monomials ea and eb. Also writes the result to ec.
function monom_lcm!(
    ec::ExponentVector{T},
    ea::ExponentVector{T},
    eb::ExponentVector{T}
) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    @inbounds ec[1] = zero(T)
    @inbounds for i in 2:length(ec)
        ec[i] = max(ea[i], eb[i])
        ec[1] += ec[i]
    end
    @invariant _monom_overflow_check(ec)
    ec
end

# Checks if the gcd of monomials ea and eb is constant.
function is_gcd_const(ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    @invariant length(ea) == length(eb)
    @inbounds for i in 2:length(ea)
        if !iszero(ea[i]) && !iszero(eb[i])
            return false
        end
    end
    true
end

# Returns the product of monomials ea and eb. Also writes the result to ec.
function monom_product!(
    ec::ExponentVector{T},
    ea::ExponentVector{T},
    eb::ExponentVector{T}
) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] + eb[j]
    end
    @invariant _monom_overflow_check(ec)
    ec
end

# Returns the result of monomial division ea / eb. Also writes the result to ec.
function monom_division!(
    ec::ExponentVector{T},
    ea::ExponentVector{T},
    eb::ExponentVector{T}
) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        @invariant ea[j] >= eb[j]
        ec[j] = ea[j] - eb[j]
    end
    ec
end

# Checks if monomial eb divides monomial ea.
function is_monom_divisible(ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    @invariant length(ea) == length(eb)
    @inbounds for j in 1:length(ea)
        if ea[j] < eb[j]
            return false
        end
    end
    true
end

# Checks if monomial eb divides monomial ea. Also writes the resulting divisor
# to ec.
function is_monom_divisible!(
    ec::ExponentVector{T},
    ea::ExponentVector{T},
    eb::ExponentVector{T}
) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        if ea[j] < eb[j]
            return false, ec
        end
        ec[j] = ea[j] - eb[j]
    end
    true, ec
end

# Checks monomials for elementwise equality
function is_monom_elementwise_eq(ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    @invariant length(ea) == length(eb)
    @inbounds for i in 1:length(ea)
        if ea[i] != eb[i]
            return false
        end
    end
    return true
end

###
# Monomial division masks.
# See f4/hashtable.jl for details.

# Constructs and returns the division mask of type Mask of the given monomial.
function monom_divmask(
    e::ExponentVector{T},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, Mask}
    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    @inbounds for i in 1:ndivvars
        for _ in 1:ndivbits
            if e[i + 1] >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end
    res
end
