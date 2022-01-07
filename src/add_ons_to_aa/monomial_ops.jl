
#------------------------------------------------------------------------------

# returns an iterator to exponent
function term_exponent_iter(m, ::Type{Val{:degrevlex}})
    N = size(m.exps, 1)
    (m.exps[j, 1] for j in 1:N - 1)
end

# returns an iterator to exponent
function term_exponent_iter(m, ::Type{Val{:lex}})
    N = size(m.exps, 1)
    (m.exps[j, 1] for j in N:-1:1)
end

#------------------------------------------------------------------------------

# returns the order of indices for traversing the monomial exponent
function exponent_ordering(R::MPolyRing, N)
    if ordering(R) == :lex
        ord = N:-1:1
    elseif ordering(R) == :degrevlex
        ord = 1:1:N - 1
    end
    ord
end

# checks if term m is divisible by the leading term of H
function is_term_divisible(m, H)
    N = size(m.exps, 1)

    for j in 1:N
        if m.exps[j, 1] < H.exps[j, 1]
            return false
        end
    end

    return true
end

# return the result of division of term m by the leading term of H,
# assumes the result is nonzero
function unsafe_term_divide(m::MPoly{T}, H::MPoly{T}) where {T}
    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    diff = MacaulayMatrix{UInt64}(undef, N, 1)
    for i in 1:N
        diff[i] = m.exps[i, 1] - H.exps[i, 1]
    end
    qu = coeff(m, 1) // coeff(H, 1)

    return MPoly{T}(R, [qu], diff)
end

#------------------------------------------------------------------------------

# return the result of dividing term m by the leading term of H
# while checking that the resulting quotiont is nonzero
function term_divides(m, H)
    R = parent(m)

    flag = is_term_divisible(m, H)
    qu = R(0)

    if flag
        qu = unsafe_term_divide(m, H)
    end

    flag, qu
end

#------------------------------------------------------------------------------

# returns gcd of leading monomials of m and H
# The coefficient in result is always set to 1
function term_gcd(m, H)
    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    supp = Int[min(m.exps[j, 1], H.exps[j, 1]) for j in ord]
    qu = R( 1 )
    set_exponent_vector!(qu, 1, supp)

    return qu
end

# checks if gcd of leading monomials of m and H is const
function is_term_gcd_constant(m, H)
    (iszero(m) || iszero(H)) && return false

    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    for j in ord
        if m.exps[j, 1] * H.exps[j, 1] > 0
            return false
        end
    end

    return true
end

#------------------------------------------------------------------------------

# returns lcm of leading monomials of m and H
# The coefficient in result is always set to 1
function term_lcm(m, H)
    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    supp = Int[max(m.exps[j, 1], H.exps[j, 1]) for j in ord]
    qu = R( 1 )
    set_exponent_vector!(qu, 1, supp)

    return qu
end

#------------------------------------------------------------------------------

# checks whether leading terms of m and H
# are equal
# Several times faster than ==
function term_equal(m, H)
    N = size(m.exps, 1)

    for j in 1:N
        if m.exps[j, 1] != H.exps[j, 1]
            return false
        end
    end

    return coeff(m, 1) == coeff(H, 1)
end

#------------------------------------------------------------------------------

# returns m*H assuming m is a single term
function multiplie_by_term(m, H::MPoly{T}) where {T}
    ans = MPoly{T}(parent(H), copy(H.coeffs), copy(H.exps))
    A = ans.exps
    for iterm in 1:size(A, 2)
        for iexp in 1:size(A, 1)
            A[iexp, iterm] += m.exps[iexp, 1]
        end
        # assuming m.coeffs[1] != 0
        ans.coeffs[iterm] *= m.coeffs[1]
    end
    ans
end
