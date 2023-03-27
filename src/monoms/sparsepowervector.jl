# SparsePowerVector{T} is a sparse vector of integers {T}
# implementing exponent vector interface.

# SparsePowerVector stores total degree in addition to partial degrees.

struct SparsePowerVector{T<:Unsigned}
    #=
    x_1^3 x_20 x_43^2  
    <-->  
    indices = [1, 20, 43]
    degrees = [3, 1, 2]
    =#
    #
    indices::Vector{UInt8}
    degrees::Vector{T}
end

# checks that there is no risk of overflow for `e`.
# If overflow if probable, throws.
function _overflow_check(e::SparsePowerVector{T}) where {T}
    _overflow_check(totaldeg(e), T)
end

# an object of type T can store capacity(T) integers at max
capacity(::Type{SparsePowerVector{T}}) where {T} = 2^32
capacity(p::SparsePowerVector{T}) where {T} = capacity(typeof(p))

# total degree of exponent vector
totaldeg(pv::PowerVector) = @inbounds pv[1]

powertype(::Type{PowerVector{T}}) where {T} = MonomHash
powertype(::PowerVector{T}) where {T} = MonomHash

make_zero_ev(::Type{PowerVector{T}}, n::Integer) where {T} = zeros(T, n+1)

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
        h += MH(x[i])*b[i]
    end
    mod(h, MonomHash)
end

#------------------------------------------------------------------------------
# Monomial orderings implementations 
# for the `PowerVector` monomial implementation.
# See monoms/orderings.jl for details.

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

# Lex exponent vector comparison
function monom_isless(ea::PowerVector, eb::PowerVector, ::Lex)
    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i] ? true : false
end

# Weighted Lex exponent vector comparison
# (Currently not exported)
function monom_isless(ea::PowerVector, eb::PowerVector, w::Weighted)
    i = 2
    weights = w.weights
    @inbounds while i < length(ea) && weights[i]*ea[i] == weights[i]*eb[i]
        i += 1
    end
    @inbounds return weights[i]*ea[i] < weights[i]*eb[i] ? true : false
end

#------------------------------------------------------------------------------
# Monomial-Monomial arithmetic.

# Returns the lcm of monomials ea and eb.
# Also writes the result to ec.
function monom_lcm!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @assert length(ec) == length(ea) == length(eb)
    @inbounds ec[1] = zero(T)
    @inbounds for i in 2:length(ec)
        ec[i] = max(ea[i], eb[i])
        ec[1] += ec[i]
    end
    ec
end

# Checks if the gcd of monomials ea and eb is constant.
function is_gcd_const(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @assert length(ea) == length(eb)
    @inbounds for i in 2:length(ea)
        if !iszero(ea[i]) && !iszero(eb[i])
            return false
        end
    end
    return true
end

# Returns the product of monomials ea and eb.
# Also writes the result to ec.
function monom_product!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @assert length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] + eb[j]
    end
    ec
end

# Returns the result of monomial division ea / eb.
# Also writes the result to ec.
function monom_division!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @assert length(ec) == length(ea) == length(eb)
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
function is_monom_divisible!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
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
    ndivvars, divmap,
    ndivbits) where {T, Mask}

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
