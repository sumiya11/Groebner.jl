# SparseExponentVector{T, N} is a sparse vector of integers of type T that
# implements the interface of a monomial of N variables.

struct SparseExponentVector{T <: Integer, N, I <: Integer}
    #=
    A monomial x_1^3 x_20 x_43^2 is stored as two vectors:
    inds = [1, 20, 43]
    vals = [3, 1, 2]
    =#
    inds::Vector{I}
    vals::Vector{T}

    function SparseExponentVector{T, N, I}(inds::Vector{I}, vals::Vector{T}) where {T, N, I}
        @invariant length(inds) == length(vals)
        new{T, N, I}(inds, vals)
    end
end

function resize_monom_if_needed!(ev, n)
    if length(ev.inds) < n
        resize!(ev.inds, n)
        resize!(ev.vals, n)
    elseif length(ev.inds) > n
        resize!(ev.inds, n)
        resize!(ev.vals, n)
    end
    nothing
end

max_vars_in_monom(::Type{SparseExponentVector{T}}) where {T} = typemax(Int8)
max_vars_in_monom(::Type{SparseExponentVector{T, N, I}}) where {T, N, I} = typemax(I)
max_vars_in_monom(p::SparseExponentVector) = max_vars_in_monom(typeof(p))

function totaldeg(ev::SparseExponentVector{T}) where {T}
    s = zero(T)
    for v in ev.vals
        s += v
    end
    s
end

copy_monom(ev::SparseExponentVector{T, N, I}) where {T, N, I} =
    SparseExponentVector{T, N, I}(Base.copy(ev.inds), Base.copy(ev.vals))

construct_const_monom(::Type{SparseExponentVector{T}}, n) where {T} =
    construct_const_monom(SparseExponentVector{T, n, Int8}, n)

function construct_const_monom(::Type{SparseExponentVector{T, N, I}}, n) where {T, N, I}
    SparseExponentVector{T, N, I}(Vector{I}(), Vector{T}())
end

construct_monom(::Type{SparseExponentVector{T}}, ev::Vector{U}) where {T, U} =
    construct_monom(SparseExponentVector{T, length(ev), Int8}, ev)

function construct_monom(
    ::Type{SparseExponentVector{T, N, I}},
    ev::Vector{U}
) where {T, N, I, U}
    inds = map(I, findall(!iszero, ev))
    vals = Vector{T}(undef, length(inds))
    @inbounds for i in 1:length(inds)
        c = ev[inds[i]]
        _monom_overflow_check(c, T)
        vals[i] = c
    end
    SparseExponentVector{T, N, I}(inds, vals)
end

function monom_to_dense_vector!(tmp::Vector{M}, pv::SparseExponentVector{T}) where {M, T}
    inds = pv.inds
    vals = pv.vals
    tmp .= zero(M)
    @inbounds for i in 1:length(inds)
        tmp[inds[i]] = vals[i]
    end
    tmp
end

construct_hash_vector(::Type{SparseExponentVector{T}}, n::Integer) where {T} =
    construct_hash_vector(SparseExponentVector{T, n, Int8}, n)

function construct_hash_vector(
    ::Type{SparseExponentVector{T, N, I}},
    n::Integer
) where {T, N, I}
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
function is_supported_ordering(
    ::Type{SparseExponentVector{T, N, I}},
    ::O
) where {T, N, I, O <: Union{Lex, DegLex, DegRevLex, InputOrdering}}
    true
end

function is_supported_ordering(
    ::Type{SparseExponentVector{T, N, I}},
    ::O
) where {T, N, I, O <: AbstractMonomialOrdering}
    false
end

# DegRevLex exponent vector comparison
function monom_isless(ea::SparseExponentVector, eb::SparseExponentVector, ::DegRevLex)
    tda, tdb = totaldeg(ea), totaldeg(eb)
    if tda < tdb
        return true
    elseif tda != tdb
        return false
    end

    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    i, j = length(ainds), length(binds)
    # TODO: @inbounds
    while i >= 1 && j >= 1
        if ainds[i] < binds[j]
            return false
        elseif ainds[i] > binds[j]
            return true
        else
            if avals[i] <= bvals[j]
                return false
            else
                return true
            end
        end
        i -= 1
        j -= 1
    end
    length(ainds) > length(binds)
end

# DegLex exponent vector comparison
function monom_isless(ea::SparseExponentVector, eb::SparseExponentVector, ::DegLex)
    tda, tdb = totaldeg(ea), totaldeg(eb)
    if tda < tdb
        return true
    elseif tda != tdb
        return false
    end
    monom_isless(ea, eb, Lex())
end

# Lex exponent vector comparison
function monom_isless(ea::SparseExponentVector, eb::SparseExponentVector, ::Lex)
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    i = 1
    while i <= length(ainds) && i <= length(binds)
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

###
# Monomial-Monomial arithmetic

function is_monom_const(ea::SparseExponentVector{T}) where {T}
    isempty(ea.inds)
end

# Returns the lcm of monomials ea and eb.
# Also writes the result to ec.
function monom_lcm!(
    ec::SparseExponentVector{T},
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    an, bn = length(ea.inds), length(eb.inds)
    resize_monom_if_needed!(ec, an + bn)
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    cinds, cvals = ec.inds, ec.vals
    @invariant length(cinds) >= an + bn
    i, j, k = 1, 1, 1
    while i <= an && j <= bn
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
    while i <= an
        cvals[k] = avals[i]
        cinds[k] = ainds[i]
        i += 1
        k += 1
    end
    while j <= bn
        cvals[k] = bvals[j]
        cinds[k] = binds[j]
        j += 1
        k += 1
    end
    resize_monom_if_needed!(ec, k - 1)
    ec
end

# Checks if the gcd of monomials ea and eb is constant.
function is_gcd_const(ea::SparseExponentVector{T}, eb::SparseExponentVector{T}) where {T}
    ainds, binds = ea.inds, eb.inds
    i, j = 1, 1
    an, bn = length(ainds), length(binds)
    # TODO: @inbounds
    while i <= an && j <= bn
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
    ec::SparseExponentVector{T},
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    # # make sure that the last index in ea is smaller than that in eb
    # if length(ea.inds) > length(eb.inds)
    #     eb, ea = ea, eb
    # end
    # @invariant length(ea.inds) <= length(eb.inds)
    an, bn = length(ea.inds), length(eb.inds)
    resize_monom_if_needed!(ec, an + bn)
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    cinds, cvals = ec.inds, ec.vals
    @invariant length(cinds) >= an + bn
    i, j, k = 1, 1, 1
    while i <= an && j <= bn
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
    while i <= an
        cvals[k] = avals[i]
        cinds[k] = ainds[i]
        i += 1
        k += 1
    end
    while j <= bn
        cvals[k] = bvals[j]
        cinds[k] = binds[j]
        j += 1
        k += 1
    end
    resize_monom_if_needed!(ec, k - 1)
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
    resize_monom_if_needed!(ec, max(an, bn))
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    cinds, cvals = ec.inds, ec.vals
    @invariant length(cinds) >= max(an, bn)
    i, j, k = 1, 1, 1
    while i <= an && j <= bn
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
    while i <= an
        cvals[k] = avals[i]
        cinds[k] = ainds[i]
        i += 1
        k += 1
    end
    resize_monom_if_needed!(ec, k - 1)
    ec
end

# Checks if monomial eb divides monomial ea.
function is_monom_divisible(
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    ainds, binds = ea.inds, eb.inds
    avals, bvals = ea.vals, eb.vals
    i, j, k = 1, 1, 1
    an, bn = length(ainds), length(binds)
    while i <= an && j <= bn
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
function is_monom_divisible!(
    ec::SparseExponentVector{T},
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    flag = is_monom_divisible(ea, eb)
    if flag
        monom_division!(ec, ea, eb)
    end
    flag, ec
end

# Checks monomials for elementwise equality
function is_monom_elementwise_eq(
    ea::SparseExponentVector{T},
    eb::SparseExponentVector{T}
) where {T}
    ea.inds == eb.inds && ea.vals == eb.vals
end

###
# Monomial division masks.
# See f4/hashtable.jl for details.

# Constructs and returns the division mask of the given monomial.
function monom_divmask(
    e::SparseExponentVector{T, N},
    DM::Type{Mask},
    ndivvars,
    divmap,
    ndivbits
) where {T, N, Mask}
    tmp = Vector{T}(undef, N)
    monom_to_dense_vector!(tmp, e)
    pushfirst!(tmp, totaldeg(e))
    monom_divmask(tmp, DM, ndivvars, divmap, ndivbits)
end
