# This file is a part of Groebner.jl. License is GNU GPL v2.

# Sparse monomial implementation
#
# SparseExponentVector{T, N} is a sparse vector of integers of type T that
# implements the interface of a monomial of N variables.

const _default_index_type = Int32

struct SparseExponentVector{T <: Integer, I <: Integer, N}
    #=
    A monomial x_1^3 x_20 x_43^2 is stored as two vectors:
    inds = [1, 20, 43]
    vals = [3, 1, 2]
    =#
    inds::Vector{I}
    vals::Vector{T}

    function SparseExponentVector{T, I, N}(inds::Vector{I}, vals::Vector{T}) where {T, I, N}
        @invariant length(inds) == length(vals)
        @invariant length(inds) <= N
        new{T, I, N}(inds, vals)
    end
end

function monom_resize_if_needed!(ev, n)
    @invariant length(ev.inds) == length(ev.vals)
    if length(ev.inds) < n
        resize!(ev.inds, n)
        resize!(ev.vals, n)
    elseif length(ev.inds) > n
        resize!(ev.inds, n)
        resize!(ev.vals, n)
    end
    nothing
end

monom_max_vars(::Type{SparseExponentVector{T}}) where {T} = typemax(_default_index_type)
monom_max_vars(::Type{SparseExponentVector{T, I, N}}) where {T, I, N} = typemax(I)
monom_max_vars(p::SparseExponentVector) = monom_max_vars(typeof(p))

function monom_totaldeg(ev::SparseExponentVector{T}) where {T}
    s = zero(T)
    for v in ev.vals
        s += v
        _monom_overflow_check(s, T)
    end
    s
end

monom_copy(ev::SparseExponentVector{T, I, N}) where {T, I, N} =
    SparseExponentVector{T, I, N}(Base.copy(ev.inds), Base.copy(ev.vals))

function monom_construct_const_monom(::Type{SparseExponentVector{T}}, n) where {T}
    monom_construct_const_monom(SparseExponentVector{T, _default_index_type, n}, n)
end

function monom_construct_const_monom(
    ::Type{SparseExponentVector{T, I, N}},
    n
) where {T, I, N}
    SparseExponentVector{T, I, N}(Vector{I}(), Vector{T}())
end

function monom_construct_from_vector(
    ::Type{SparseExponentVector{T}},
    ev::Vector{U}
) where {T, U}
    n = length(ev)
    monom_construct_from_vector(SparseExponentVector{T, _default_index_type, n}, ev)
end

function monom_construct_from_vector(
    ::Type{SparseExponentVector{T, I, N}},
    ev::Vector{U}
) where {T, I, N, U}
    @invariant length(ev) == N
    inds = map(I, findall(!iszero, ev))
    vals = Vector{T}(undef, length(inds))
    @inbounds for i in 1:length(inds)
        c = ev[inds[i]]
        _monom_overflow_check(c, T)
        vals[i] = c
    end
    SparseExponentVector{T, I, N}(inds, vals)
end

function monom_to_vector!(tmp::Vector{M}, pv::SparseExponentVector{T}) where {M, T}
    inds = pv.inds
    vals = pv.vals
    tmp .= zero(M)
    @inbounds for i in 1:length(inds)
        tmp[inds[i]] = vals[i]
    end
    tmp
end

function monom_construct_hash_vector(::Type{SparseExponentVector{T}}, n::Integer) where {T}
    monom_construct_hash_vector(SparseExponentVector{T, _default_index_type, n}, n)
end

function monom_construct_hash_vector(
    ::Type{SparseExponentVector{T, I, N}},
    n::Integer
) where {T, I, N}
    @invariant n == N
    rand(MonomHash, n)
end

function monom_hash(ev::SparseExponentVector{T}, b::Vector{MH}) where {T, MH}
    h = zero(MH)
    inds = ev.inds
    vals = ev.vals
    @inbounds for i in 1:length(inds)
        idx = inds[i]
        h += MH(vals[i]) * b[idx]
    end
    mod(h, MonomHash)
end

###
# Monomial comparator functions. See monoms/orderings.jl for details.

# SparseExponentVector supports only lex, deglex, and degrevlex
function monom_is_supported_ordering(
    ::Type{SparseExponentVector{T, I, N}},
    ::O
) where {T, I, N, O}
    O <: Union{Lex, DegLex, DegRevLex, InputOrdering}
end

# DegRevLex exponent vector comparison
function monom_isless(
    ea::SparseExponentVector,
    eb::SparseExponentVector,
    ::_DegRevLex{true}
)
    tda, tdb = monom_totaldeg(ea), monom_totaldeg(eb)
    if tda < tdb
        return true
    elseif tda != tdb
        return false
    end

    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    i, j = length(ainds), length(binds)
    @inbounds while i >= 1 && j >= 1
        if ainds[i] < binds[j]
            return false
        elseif ainds[i] > binds[j]
            return true
        else
            if avals[i] < bvals[j]
                return false
            elseif avals[i] > bvals[j]
                return true
            end
        end
        i -= 1
        j -= 1
    end
    length(ainds) > length(binds)
end

# DegLex exponent vector comparison
function monom_isless(ea::SparseExponentVector, eb::SparseExponentVector, ::_DegLex{true})
    tda, tdb = monom_totaldeg(ea), monom_totaldeg(eb)
    if tda < tdb
        return true
    elseif tda != tdb
        return false
    end
    monom_isless(ea, eb, _Lex{true}(Int[]))
end

# Lex exponent vector comparison
function monom_isless(ea::SparseExponentVector, eb::SparseExponentVector, ::_Lex{true})
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    i = 1
    @inbounds while i <= length(ainds) && i <= length(binds)
        if ainds[i] < binds[i]
            return false
        elseif ainds[i] > binds[i]
            return true
        else
            if avals[i] < bvals[i]
                return true
            elseif avals[i] > bvals[i]
                return false
            end
        end
        i += 1
    end
    length(ainds) < length(binds)
end

function monom_isless(
    ea::SparseExponentVector{T, I, N},
    eb::SparseExponentVector{T, I, N},
    ord
) where {T, I, N}
    tmp1, tmp2 = Vector{T}(undef, N), Vector{T}(undef, N)
    a = monom_construct_from_vector(ExponentVector{T}, monom_to_vector!(tmp1, ea))
    b = monom_construct_from_vector(ExponentVector{T}, monom_to_vector!(tmp2, eb))
    monom_isless(a, b, ord)
end

###
# Monomial-Monomial arithmetic

function is_monom_const(ea::SparseExponentVector{T}) where {T}
    isempty(ea.inds)
end

# Returns the lcm of monomials ea and eb.
# Also writes the result to ec.
function monom_lcm!(
    ec::SparseExponentVector{T, I, N},
    ea::SparseExponentVector{T, I, N},
    eb::SparseExponentVector{T, I, N}
) where {T, I, N}
    an, bn = length(ea.inds), length(eb.inds)
    monom_resize_if_needed!(ec, an + bn)
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    cinds, cvals = ec.inds, ec.vals
    @invariant length(cinds) >= an + bn
    i, j, k = 1, 1, 1
    @inbounds while i <= an && j <= bn
        ai, bj = ainds[i], binds[j]
        if ai == bj
            cvals[k] = max(avals[i], bvals[j])
            cinds[k] = ai
            i += 1
            j += 1
        elseif ai < bj
            cvals[k] = avals[i]
            cinds[k] = ai
            i += 1
        else
            cvals[k] = bvals[j]
            cinds[k] = bj
            j += 1
        end
        k += 1
    end
    @inbounds while i <= an
        cvals[k] = avals[i]
        cinds[k] = ainds[i]
        i += 1
        k += 1
    end
    @inbounds while j <= bn
        cvals[k] = bvals[j]
        cinds[k] = binds[j]
        j += 1
        k += 1
    end
    monom_resize_if_needed!(ec, k - 1)
    ec
end

# Checks if the gcd of monomials ea and eb is constant.
function monom_is_gcd_const(
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    ainds, binds = ea.inds, eb.inds
    i, j = 1, 1
    an, bn = length(ainds), length(binds)
    @inbounds while i <= an && j <= bn
        if ainds[i] == binds[j]
            return false
        end
        if ainds[i] < binds[j]
            i += 1
        else
            j += 1
        end
    end
    true
end

# Returns the product of monomials ea and eb.
# Also writes the result to ec.
function monom_product!(
    ec::SparseExponentVector{T, I, N},
    ea::SparseExponentVector{T, I, N},
    eb::SparseExponentVector{T, I, N}
) where {T, I, N}
    an, bn = length(ea.inds), length(eb.inds)
    monom_resize_if_needed!(ec, min(N, an + bn))

    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    cinds, cvals = ec.inds, ec.vals

    i, j, k = 1, 1, 1
    @inbounds while i <= an && j <= bn
        ai, bj = ainds[i], binds[j]
        if ai == bj
            cvals[k] = avals[i] + bvals[j]
            _monom_overflow_check(cvals[k], T)
            @invariant !iszero(cvals[k])
            cinds[k] = ai
            i += 1
            j += 1
        elseif ai < bj
            cvals[k] = avals[i]
            cinds[k] = ai
            i += 1
        else
            cvals[k] = bvals[j]
            cinds[k] = bj
            j += 1
        end
        k += 1
    end
    @inbounds while i <= an
        cvals[k] = avals[i]
        cinds[k] = ainds[i]
        i += 1
        k += 1
    end
    @inbounds while j <= bn
        cvals[k] = bvals[j]
        cinds[k] = binds[j]
        j += 1
        k += 1
    end

    monom_resize_if_needed!(ec, k - 1)
    ec
end

# Returns the result of monomial division ea / eb.
# Also writes the result to ec.
function monom_division!(
    ec::SparseExponentVector{T},
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    an, bn = length(ea.inds), length(eb.inds)
    monom_resize_if_needed!(ec, max(an, bn))
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    cinds, cvals = ec.inds, ec.vals
    @invariant length(cinds) >= max(an, bn)
    i, j, k = 1, 1, 1
    @inbounds while i <= an && j <= bn
        ai, bj = ainds[i], binds[j]
        if ai == bj
            @invariant avals[i] >= bvals[j]
            d = avals[i] - bvals[j]
            if !iszero(d)
                cvals[k] = d
                cinds[k] = ai
                k += 1
            end
            i += 1
            j += 1
        else
            @invariant ai < bj
            cvals[k] = avals[i]
            cinds[k] = ai
            i += 1
            k += 1
        end
    end
    @invariant j - 1 == bn
    @inbounds while i <= an
        cvals[k] = avals[i]
        cinds[k] = ainds[i]
        i += 1
        k += 1
    end
    monom_resize_if_needed!(ec, k - 1)
    ec
end

# Checks if monomial eb divides monomial ea.
function monom_is_divisible(
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    i, j, k = 1, 1, 1
    an, bn = length(ainds), length(binds)
    @inbounds while i <= an && j <= bn
        ai, bj = ainds[i], binds[j]
        if ai == bj
            if avals[i] < bvals[j]
                return false
            end
            i += 1
            j += 1
        elseif ai < bj
            i += 1
        else
            return false
        end
        k += 1
    end
    i <= an || (i - 1 == an && j - 1 == bn)
end

# Checks if monomial eb divides monomial ea.
# Also writes the resulting divisor to ec.
function monom_is_divisible!(
    ec::SparseExponentVector{T},
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    flag = monom_is_divisible(ea, eb)
    if flag
        monom_division!(ec, ea, eb)
    end
    flag, ec
end

# Checks monomials for elementwise equality
function monom_is_equal(ea::SparseExponentVector{T}, eb::SparseExponentVector{T}) where {T}
    length(ea.inds) != length(eb.inds) && return false
    @inbounds for i in 1:length(ea.inds)
        !(ea.inds[i] == eb.inds[i] && ea.vals[i] == eb.vals[i]) && return false
    end
    true
end

###
# Monomial division masks.
# See f4/hashtable.jl for details.

# Constructs and returns the division mask of the given monomial.
function monom_create_divmask(
    e::SparseExponentVector{T, I, N},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, I, N, Mask}
    tmp = Vector{T}(undef, N)
    monom_to_vector!(tmp, e)
    pushfirst!(tmp, monom_totaldeg(e))
    monom_create_divmask(tmp, DM, ndivvars, divmap, ndivbits)
end
