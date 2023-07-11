
# PackedVector{T, B} admits the interface of a vector
# of small nonegative integers with O(1) vector sum. 
# 
# PackedVector{T, B} is implemented as a list of chunks,
# each chunk is an integer {T} that packs a fixed number
# of small numbers {B} together.
#
# PackedVector always stores the sum of vector explicitly.
# PackedVector always stores in degree-reverse-lex favorable order.
#
# if nchunks = 2 and chunksize = 3, then the sequence
# of numbers (3, 2, 7, 6) is stored as a pair of chunks:
# (18 * 2^16 + 6 * 2^8 + 7, 2 * 2^16 + 3 * 2^8 + 0)
#
# conversion from dense vector to PackedVector{T, B}: free
# conversion from PackedVector{T, B} to dense vector: need length
struct PackedVector{T <: Unsigned, B <: Unsigned}
    data::Vector{T}
end

@noinline function _monom_overflow_check(e::PackedVector{T, B}) where {T, B}
    totaldeg(e) >= div(typemax(B), 2) && throw("Overflow is probable.")
end

Base.eltype(::Type{PackedVector{T, B}}) where {T, B} = B
Base.eltype(::PackedVector{T, B}) where {T, B} = B
Base.copy(pv::PackedVector{T, B}) where {T, B} = PackedVector{T, B}(copy(pv.data))

max_vars_in_monom(::Type{PackedVector{T, B}}) where {T, B} = 2^32
max_vars_in_monom(::PackedVector{T, B}) where {T, B} = max_vars_in_monom(PackedVector{T, B})

PackedVector(ev) = PackedVector{UInt64}(ev)
PackedVector{T}(ev) where {T <: Unsigned} = PackedVector{T, UInt8}(ev)

function totaldeg(pv::PackedVector{T, B}) where {T, B}
    B(@inbounds(pv.data[1]) >> (8 * (sizeof(T) - sizeof(B))))
end

function construct_monom(::Type{PackedVector{T, B}}, ev::Vector{U}) where {T, B, U}
    ts, bs = sizeof(T), sizeof(B)
    epc = div(ts, bs)
    @assert epc * bs == ts
    n = length(ev)
    nch = nchunks(T, B, n)
    indent = ts - degsize(T, B, n)
    data = zeros(T, nch)
    s = zero(T)
    @inbounds for i in n:-1:1
        d = T(ev[i])
        j = nch + 1 - (div(i - 1, epc) + 1)
        islast = iszero(mod(i - 1, epc))
        data[j] = data[j] | d
        !islast && (data[j] = data[j] << (bs * 8))
        s += d
    end
    @inbounds data[1] |= s << (indent * 8)
    PackedVector{T, B}(data)
end

function make_zero_ev(::Type{PackedVector{T, B}}, n::Integer) where {T, B}
    PackedVector{T, B}(zeros(T, nchunks(T, B, n)))
end

function make_dense(pv::PackedVector{T, B}, n::Integer) where {T, B}
    # @warn "Make dense is called"
    d = Vector{B}(undef, n)
    epc = div(sizeof(T), sizeof(B))
    k = 0
    for c in length(pv.data):-1:1
        idx1 = k + 1
        idx2 = min(k + epc, n)
        dd = idx2 - idx1 + 1
        d[idx1:idx2] .= decompose(pv.data[c], B)[1:dd]
        k += epc
    end
    d
end

function make_hasher(::Type{PackedVector{T, B}}, n::Integer) where {T, B}
    PackedVector{T, B}(rand(T, nchunks(T, B, n)))
end

function Base.hash(x::PackedVector{T, B}, b::PackedVector{T, B}) where {T, B}
    @assert length(x.data) == length(b.data)
    maxT = promote_type(T, MonomHash)
    h = zero(maxT)
    @inbounds for i in 1:length(x.data)
        h += maxT(x.data[i]) * b.data[i]
    end
    mod(h, MonomHash)
end

function revealify(pv::PackedVector{T, B}) where {T, B}
    ts, bs = sizeof(T), sizeof(B)
    epc = div(ts, bs)
    ans = Vector{Vector{String}}()
    for chunk in pv.data
        bs = bitstring(chunk)
        s = [bs[((i - 1) * 8 + 1):(i * 8)] for i in 1:epc]
        push!(ans, s)
    end
    ans
end

#------------------------------------------------------------------------------

function exponent_isless_lex(ea::PackedVector{T, B}, eb::PackedVector{T, B}) where {T, B}
    s = div(sizeof(T), sizeof(B)) * length(ea.data)
    a = construct_monom(ExponentVector{T}, make_dense(ea, s))
    b = construct_monom(ExponentVector{T}, make_dense(eb, s))
    exponent_isless_lex(a, b)
end

function exponent_isless_dl(ea::PackedVector{T, B}, eb::PackedVector{T, B}) where {T, B}
    s = div(sizeof(T), sizeof(B)) * length(ea.data)
    a = construct_monom(ExponentVector{T}, make_dense(ea, s))
    b = construct_monom(ExponentVector{T}, make_dense(eb, s))
    exponent_isless_dl(a, b)
end

function exponent_isless_drl(ea::PackedVector{T, B}, eb::PackedVector{T, B}) where {T, B}
    if totaldeg(ea) < totaldeg(eb)
        return true
    end
    if totaldeg(ea) > totaldeg(eb)
        return false
    end
    i = 1
    @inbounds while i < length(ea.data) && ea.data[i] == eb.data[i]
        i += 1
    end
    if ea.data[i] <= eb.data[i]
        return false
    else
        return true
    end
end

#------------------------------------------------------------------------------

function monom_lcm!(
    ec::PackedVector{T, B},
    ea::PackedVector{T, B},
    eb::PackedVector{T, B}
) where {T, B}
    @inbounds ec.data[1] = zero(T)
    @inbounds for i in 1:length(ec.data)
        a, b = ea.data[i], eb.data[i]
        ec.data[i], si = packedmax(a, b, B, Val(isone(i)))
        ec.data[1] += si << ((sizeof(T) - sizeof(B)) * 8)
    end
    _monom_overflow_check(ec)
    ec
end

function is_gcd_const(ea::PackedVector{T, B}, eb::PackedVector{T, B}) where {T, B}
    @inbounds for i in 1:length(ea.data)
        if !packedorth(ea.data[i], eb.data[i], B, Val(isone(i)))
            return false
        end
    end
    return true
end

function monom_product!(
    ec::PackedVector{T, B},
    ea::PackedVector{T, B},
    eb::PackedVector{T, B}
) where {T, B}
    @assert length(ec.data) == length(ea.data) == length(eb.data)
    for i in 1:length(ec.data)
        ec.data[i] = ea.data[i] + eb.data[i]
    end
    _monom_overflow_check(ec)
    ec
end

function monom_division!(
    ec::PackedVector{T, B},
    ea::PackedVector{T, B},
    eb::PackedVector{T, B}
) where {T, B}
    @assert length(ec.data) == length(ea.data) == length(eb.data)
    for i in 1:length(ec.data)
        ec.data[i] = ea.data[i] - eb.data[i]
    end
    ec
end

function is_monom_divisible(ea::PackedVector{T, B}, eb::PackedVector{T, B}) where {T, B}
    @inbounds for i in 1:length(ea.data)
        if !packedge(ea.data[i], eb.data[i], B, Val(isone(i)))
            return false
        end
    end
    return true
end

function is_monom_divisible!(
    ec::PackedVector{T, B},
    ea::PackedVector{T, B},
    eb::PackedVector{T, B}
) where {T, B}
    ans = is_monom_divisible(ea, eb)
    ans && monom_division!(ec, ea, eb)
    ans, ec
end

function is_monom_elementwise_eq(
    ea::PackedVector{T, B},
    eb::PackedVector{T, B}
) where {T, B}
    @inbounds ea.data[1] != eb.data[1] && return false
    @inbounds for i in 2:length(ea.data)
        if ea.data[i] != eb.data[i]
            return false
        end
    end
    return true
end

#------------------------------------------------------------------------------

function monom_divmask(
    e::PackedVector{T, B},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, B, Mask}
    ee = make_dense(e, length(e.data) * div(sizeof(T), sizeof(B)))

    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    for i in 1:ndivvars
        for j in 1:ndivbits
            @inbounds if ee[i] >= divmap[ctr]
                res |= o << (ctr - 1)
            end
            ctr += o
        end
    end

    res
end
