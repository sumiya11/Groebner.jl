# PowerVector{T} is a dense vector of integers {T}
# implementing exponent vector interface.

# PowerVector stores total degree in addition to partial degrees.

const PowerVector{T} = Vector{T} where {T<:Integer}

# checks that there is no risk of overflow for `e`.
# If overflow if probable, throws.
function _overflow_check(e::PowerVector{T}) where {T}
    _overflow_check(totaldeg(e), T)
end

# an object of type T can store capacity(T) integers at max
capacity(::Type{PowerVector{T}}) where {T} = 2^32
capacity(p::PowerVector{T}) where {T} = capacity(typeof(p))

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

#=
    :degrevlex exponent vector comparison
=#
function exponent_isless_drl(ea::PowerVector, eb::PowerVector)
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

    if ea[i] <= eb[i]
        return false
    else
        return true
    end
end

#=
    :deglex exponent vector comparison
=#
function exponent_isless_dl(ea::PowerVector, eb::PowerVector)
    if @inbounds ea[1] < eb[1]
        return true
    elseif @inbounds ea[1] != eb[1]
        return false
    end

    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    return ea[i] < eb[i] ? true : false
end

#=
    :lex exponent vector comparison
=#
function exponent_isless_lex(ea::PowerVector, eb::PowerVector)
    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    return ea[i] < eb[i] ? true : false
end

#------------------------------------------------------------------------------

function monom_lcm!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds ec[1] = zero(T)
    @inbounds for i in 2:length(ec)
        ec[i] = max(ea[i], eb[i])
        ec[1] += ec[i]
    end
    ec
end

function is_gcd_const(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for i in 2:length(ea)
        if !iszero(ea[i]) && !iszero(eb[i])
            return false
        end
    end
    return true
end

function monom_product!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] + eb[j]
    end
    ec
end

function monom_division!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] - eb[j]
    end
    ec
end

function is_monom_divisible(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for j in 1:length(ea)
        if ea[j] < eb[j]
            return false
        end
    end
    return true
end

function is_monom_divisible!(ec::PowerVector{T}, ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for j in 1:length(ec)
        if ea[j] < eb[j]
            return false, ec
        end
        ec[j] = ea[j] - eb[j]
    end
    return true, ec
end

function is_monom_elementwise_eq(ea::PowerVector{T}, eb::PowerVector{T}) where {T}
    @inbounds for i in 1:length(ea)
        if ea[i] != eb[i]
            return false
        end
    end
    return true
end

#------------------------------------------------------------------------------

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
