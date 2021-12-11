
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
    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    for j in ord
        if m.exps[j, 1] < H.exps[j, 1]
            return false
        end
    end

    return true
end

# return the result of division of term m by the leading term of H,
# assumes the result is nonzero
function unsafe_term_divide(m, H)
    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    diff = Int[m.exps[j, 1] - H.exps[j, 1] for j in ord]
    qu = R( coeff(m, 1) // coeff(H, 1) )
    set_exponent_vector!(qu, 1, diff)

    return qu
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
function term_equal(m, H)
    R = parent(m)
    N = size(m.exps, 1)
    ord = exponent_ordering(R, N)

    for j in ord
        if m.exps[j, 1] != H.exps[j, 1]
            return false
        end
    end

    return coeff(m, 1) == coeff(H, 1)
end
