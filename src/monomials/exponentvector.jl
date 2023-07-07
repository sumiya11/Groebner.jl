# PowerVector{T} is a dense vector of integers {T}
# implementing exponent vector interface.

# PowerVector stores the total degree in addition to partial degrees.
# For example, x^2yz^5 may be stored as a *dynamic* vector [8, 2, 1, 5].

# PowerVector.
# PowerVector is just a dynamic vector of integers.
const PowerVector{T} = Vector{T} where {T <: Integer}

# Checks whether PowerVector{T} provides a comparator function implementation
# for the given monomial ordering of type `O`
function is_supported_ordering(::Type{PowerVector{T}}, ::O) where {T, O}
    O <: AbstractMonomialOrdering
end

# Checks that there is no risk of exponent overflow for `e`.
# If overflow if probable, throws a RecoverableException.
function _overflow_check(e::PowerVector{T}) where {T}
    _overflow_check(totaldeg(e), T)
end

# The maximum number of variables that a PowerVector can store.
capacity(::Type{PowerVector{T}}) where {T} = 2^32
capacity(p::PowerVector{T}) where {T} = capacity(typeof(p))

# The total degree of an exponent vector
totaldeg(pv::PowerVector) = @inbounds pv[1]

# The type of an entry of a PowerVector{T}. 
# Note that this is not necessarily equal to T.
# Change to entrytype?
powertype(::Type{PowerVector{T}}) where {T} = MonomHash
powertype(::PowerVector{T}) where {T} = MonomHash

make_zero_ev(::Type{PowerVector{T}}, n::Integer) where {T} = zeros(T, n + 1)

function make_ev(::Type{PowerVector{T}}, ev::Vector{U}) where {T, U}
    v = Vector{T}(undef, length(ev) + 1)
    @inbounds v[1] = sum(ev)
    @inbounds v[2:end] .= ev
    PowerVector{T}(v)
end

function make_dense!(tmp::Vector{M}, pv::PowerVector{T}) where {M, T}
    @assert length(tmp) == length(pv) - 1
    @inbounds tmp[1:end] = pv[2:end]
    tmp
end

function make_hasher(::Type{PowerVector{T}}, n::Integer) where {T}
    rand(MonomHash, n + 1)
end

function Base.hash(x::PowerVector{T}, b::PowerVector{MH}) where {T, MH}
    h = zero(MH)
    @inbounds for i in eachindex(x, b)
        h += MH(x[i]) * b[i]
    end
    mod(h, MonomHash)
end

#------------------------------------------------------------------------------
# Monomial orderings comparator functions implementations
# for the `PowerVector` monomial implementation.
# See monoms/orderings.jl for details.

#####
# DegRevLex exponent vector comparison
function monom_isless(ea::PowerVector, eb::PowerVector, ::DegRevLex)
    if @inbounds ea[1] < eb[1]
        return true
    elseif @inbounds ea[1] != eb[1]
        return false
    end
    # assuming length(ea) > 1
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
# DegRevLex, but on a given range of variables [lo:hi].
# The flag hasdegree indicates whether the total degree is stored 
# in ea[lo] and eb[lo], respectively.
function monom_isless(
    ea::PowerVector,
    eb::PowerVector,
    ::DegRevLex,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
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

#####
# DegLex exponent vector comparison
function monom_isless(ea::PowerVector, eb::PowerVector, ::DegLex)
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
# DegLex, but on a given range of variables [lo:hi].
# The flag hasdegree indicates whether the total degree is stored 
# in ea[lo] and eb[lo], respectively.
function monom_isless(
    ea::PowerVector,
    eb::PowerVector,
    ::DegLex,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
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

#####
# Lex
function monom_isless(ea::PowerVector, eb::PowerVector, ::Lex)
    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end
# Lex, but on a given range of variables [lo:hi].
# The flag hasdegree indicates whether the total degree is stored 
# in ea[lo] and eb[lo], respectively.
function monom_isless(
    ea::PowerVector,
    eb::PowerVector,
    ::Lex,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    i = lo + hasdegree
    @inbounds while i < hi && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end

#####
# Weighted exponent vector comparison.
function monom_isless(ea::PowerVector, eb::PowerVector, w::WeightedOrdering)
    monom_isless(ea, eb, w, 1, length(ea), true)
end
# Weighted, but on a given range of variables [lo:hi].
function monom_isless(
    ea::PowerVector,
    eb::PowerVector,
    w::WeightedOrdering,
    lo::Int,
    hi::Int,
    hasdegree::Bool
)
    weights = w.weights
    common_type = promote_type(eltype(weights), eltype(ea))
    sa, sb = zero(common_type), zero(common_type)
    j = 1
    @inbounds for i in (lo + hasdegree):hi
        sa += weights[j] * ea[i]
        sb += weights[j] * eb[i]
        j += 1
    end
    if sa < sb
        true
    elseif sa == sb
        monom_isless(ea, eb, w.ord, lo + hasdegree, hi, false)
    else
        false
    end
end

#####
# Block exponent vector comparison.
function monom_isless(ea::PowerVector, eb::PowerVector, b::BlockOrdering)
    monom_isless(ea, eb, b, 1, length(ea), true)
end
# Block, but on a given range of variables [lo:hi].
function monom_isless(
    ea::PowerVector,
    eb::PowerVector,
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
function monom_isless(ea::PowerVector, eb::PowerVector, m::MatrixOrdering)
    monom_isless(ea, eb, m, 1, length(ea), true)
end
# Matrix, but on a given range of variables [lo:hi].
function monom_isless(
    ea::PowerVector,
    eb::PowerVector,
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

#------------------------------------------------------------------------------
# Monomial-Monomial arithmetic.

# Returns the lcm of monomials ea and eb.
# Also writes the result to ec.
function monom_lcm!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    # @assert length(ec) == length(ea) == length(eb)
    @inbounds ec[1] = zero(T)
    @inbounds for i in 2:length(ec)
        ec[i] = max(ea[i], eb[i])
        ec[1] += ec[i]
    end
    ec
end

# Checks if the gcd of monomials ea and eb is constant.
function is_gcd_const(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    # @assert length(ea) == length(eb)
    @inbounds for i in 2:length(ea)
        if !iszero(ea[i]) && !iszero(eb[i])
            return false
        end
    end
    return true
end

# Returns the product of monomials ea and eb.
# Also writes the result to ec.
function monom_product!(
    ec::PowerVector{T},
    ea::PowerVector{T},
    eb::PowerVector{T}
) where {T}
    # @assert length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] + eb[j]
    end
    ec
end

# Returns the result of monomial division ea / eb.
# Also writes the result to ec.
function monom_division!(
    ec::PowerVector{T},
    ea::PowerVector{T},
    eb::PowerVector{T}
) where {T}
    # @assert length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] - eb[j]
    end
    ec
end

# Checks if monomial eb divides monomial ea.
function is_monom_divisible(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for j in 1:length(ea)
        if ea[j] < eb[j]
            return false
        end
    end
    return true
end

# Checks if monomial eb divides monomial ea.
# Also writes the resulting divisor to ec.
function is_monom_divisible!(
    ec::PowerVector{T},
    ea::PowerVector{T},
    eb::PowerVector{T}
) where {T}
    @inbounds for j in 1:length(ec)
        if ea[j] < eb[j]
            return false, ec
        end
        ec[j] = ea[j] - eb[j]
    end
    return true, ec
end

# Checks monomials for elementwise equality
function is_monom_elementwise_eq(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for i in 1:length(ea)
        if ea[i] != eb[i]
            return false
        end
    end
    return true
end

#------------------------------------------------------------------------------
# Monomial division masks.
# See f4/hashtable.jl for details.

# Constructs and returns the division mask of the given monomial.
function monom_divmask(
    e::PowerVector{T},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, Mask}
    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    for i in 1:ndivvars
        for j in 1:ndivbits
            @inbounds if e[i + 1] >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end
