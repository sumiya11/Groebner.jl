# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# ExponentVector implements the monomial interface (whatever this interface is).
#
# ExponentVector{T} is a vector of integers of type T. It is commonly used with
# the UInt8 element type (that is 8 bits per exponent).
#
# Note: Some functions that implement monomial arithmetic and monomial
# comparisons are not automatically vectorized with the use of AVX2 or AVX512
# (in Julia 1.10). When vectorized, the code does not handle scalar tails well.
# At the same time, in hand-vectorized code it is hard to do masked load.  

const ExponentVector{T} = Vector{T} where {T <: Integer}

# MonomialDegreeOverflow is thrown if there is a risk of monomial degree
# overflow. If we catch a MonomialDegreeOverflow, there is some hope to recover
# the program by restarting with a wider integer type for storing exponents.
struct MonomialDegreeOverflow <: Exception
    msg::String
end

Base.showerror(io::IO, e::MonomialDegreeOverflow) = print(io, e.msg)

@noinline __throw_monom_overflow_error(c, B) =
    throw(MonomialDegreeOverflow("Overflow may happen with the entry $c of type $B."))

monom_overflow_threshold(::Type{T}) where {T <: Integer} = div(typemax(T), 2)
monom_overflow_check(a::T) where {T <: Integer} = monom_overflow_check(a, T)

function monom_overflow_check(a::Integer, ::Type{T}) where {T}
    a >= monom_overflow_threshold(T) && __throw_monom_overflow_error(a, T)
    true
end

function monom_overflow_check(e::ExponentVector{T}) where {T}
    # Exponent vector overflows if the total degree overflows
    monom_overflow_check(monom_totaldeg(e))
end

monom_max_vars(::Type{ExponentVector{T}}) where {T} = 2^32
monom_max_vars(p::ExponentVector{T}) where {T} = monom_max_vars(typeof(p))

monom_totaldeg(pv::ExponentVector) = @inbounds pv[1]
monom_entrytype(pv::ExponentVector{T}) where {T} = T
monom_entrytype(::Type{ExponentVector{T}}) where {T} = T

monom_copy(pv::ExponentVector) = Base.copy(pv)

function monom_copy!(dst::ExponentVector, src::ExponentVector)
    @invariant length(dst) == length(src)
    @inbounds for i in 1:length(src)
        dst[i] = src[i]
    end
    dst
end

monom_construct_const(::Type{ExponentVector{T}}, n::Integer) where {T} = zeros(T, n + 1)

function monom_construct_from_vector(::Type{ExponentVector{T}}, ev::Vector{U}) where {T, U}
    v = ExponentVector{T}(undef, length(ev) + 1)
    s = zero(T)
    @inbounds for i in 1:length(ev)
        monom_overflow_check(ev[i], T)
        s += T(ev[i])
        monom_overflow_check(s, T)
        v[i + 1] = T(ev[i])
    end
    @inbounds v[1] = s
    v
end

function monom_to_vector!(tmp::Vector{M}, pv::ExponentVector{T}) where {M, T}
    @invariant length(tmp) == length(pv) - 1
    @inbounds tmp[1:end] = pv[2:end]
    tmp
end

function monom_construct_hash_vector(
    rng::AbstractRNG,
    ::Type{ExponentVector{T}},
    n::Integer
) where {T}
    rand(rng, MonomHash, n + 1)
end

function monom_hash(x::ExponentVector{T}, b::Vector{MH}) where {T, MH}
    @invariant length(x) == length(b)
    h = zero(MH)
    @inbounds for i in 1:length(x)
        h += MH(x[i]) * b[i]
    end
    mod(h, MonomHash)
end

###
# Monomial comparators.

function monom_is_supported_ordering(::Type{ExponentVector{T}}, ::Any) where {T}
    # ExponentVector implements all monomial orderings.
    true
end

# DegRevLex monomial comparison. 
function monom_isless(ea::ExponentVector, eb::ExponentVector, ::DegRevLex{true})
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    if monom_totaldeg(ea) < monom_totaldeg(eb)
        return true
    elseif monom_totaldeg(ea) != monom_totaldeg(eb)
        return false
    end
    i = length(ea)
    @inbounds while i > 2 && ea[i] == eb[i]
        i -= 1
    end
    @inbounds return ea[i] > eb[i]
end

# DegRevLex monomial comparison (shuffled variables).
function monom_isless(ea::ExponentVector{T}, eb::ExponentVector{T}, ord::DegRevLex{false}) where {T}
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    variables = ordering_variables(ord)
    eatotaldeg = zero(T)
    ebtotaldeg = zero(T)
    @inbounds for i in variables
        eatotaldeg += ea[i + 1]
        ebtotaldeg += eb[i + 1]
    end
    if eatotaldeg < ebtotaldeg
        return true
    elseif eatotaldeg != ebtotaldeg
        return false
    end
    i = length(variables)
    @inbounds while i > 1 && ea[variables[i] + 1] == eb[variables[i] + 1]
        i -= 1
    end
    @inbounds ea[variables[i] + 1] > eb[variables[i] + 1]
end

# DegLex monomial comparison.
function monom_isless(ea::ExponentVector, eb::ExponentVector, ::DegLex{true})
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    if monom_totaldeg(ea) < monom_totaldeg(eb)
        return true
    elseif monom_totaldeg(ea) != monom_totaldeg(eb)
        return false
    end
    monom_isless(ea, eb, Lex())
end

# DegLex monomial comparison (shuffled variables).
function monom_isless(ea::ExponentVector{T}, eb::ExponentVector{T}, ord::DegLex{false}) where {T}
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    variables = ordering_variables(ord)
    eatotaldeg = zero(T)
    ebtotaldeg = zero(T)
    @inbounds for i in variables
        eatotaldeg += ea[i + 1]
        ebtotaldeg += eb[i + 1]
    end
    if eatotaldeg < ebtotaldeg
        return true
    elseif eatotaldeg != ebtotaldeg
        return false
    end
    i = 1
    @inbounds while i < length(variables) && ea[variables[i] + 1] == eb[variables[i] + 1]
        i += 1
    end
    @inbounds ea[variables[i] + 1] < eb[variables[i] + 1]
end

# Lex monomial comparison.
function monom_isless(ea::ExponentVector, eb::ExponentVector, ::Lex{true})
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    i = 2
    @inbounds while i < length(ea) && ea[i] == eb[i]
        i += 1
    end
    @inbounds return ea[i] < eb[i]
end

# Lex monomial comparison (shuffled variables).
function monom_isless(ea::ExponentVector, eb::ExponentVector, ord::Lex{false})
    variables = ordering_variables(ord)
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    i = 1
    @inbounds while i < length(variables) && ea[variables[i] + 1] == eb[variables[i] + 1]
        i += 1
    end
    @inbounds ea[variables[i] + 1] < eb[variables[i] + 1]
end

# Weighted monomial comparison.
function monom_isless(
    ea::ExponentVector{T},
    eb::ExponentVector{T},
    ord::WeightedOrdering{U}
) where {U, T}
    weights = ord.weights
    variables = ordering_variables(ord)
    @invariant length(ea) == length(eb)
    @invariant length(ea) > 1
    sa, sb = zero(U), zero(U)
    @inbounds for i in 1:length(variables)
        sa += ea[variables[i] + 1] * weights[i]
        sb += eb[variables[i] + 1] * weights[i]
    end
    if sa < sb
        return true
    elseif sa != sb
        return false
    end
    monom_isless(ea, eb, Lex())
end

# Product ordering exponent vector comparison.
function monom_isless(ea::ExponentVector, eb::ExponentVector, b::ProductOrdering)
    if monom_isless(ea, eb, b.ord1)
        return true
    end
    if monom_isless(eb, ea, b.ord1)
        return false
    end
    monom_isless(ea, eb, b.ord2)
end

# Matrix ordering exponent vector comparison.
function monom_isless(ea::ExponentVector, eb::ExponentVector, m::MatrixOrdering)
    rows = m.rows
    variables = ordering_variables(m)
    @inbounds common_type = signed(promote_type(eltype(rows[1]), eltype(ea)))
    @inbounds for i in 1:length(rows)
        sa, sb = zero(common_type), zero(common_type)
        for j in 1:length(variables)
            sa += rows[i][j] * common_type(ea[variables[j] + 1])
            sb += rows[i][j] * common_type(eb[variables[j] + 1])
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
# Monomial-monomial arithmetic.

# Returns the lcm of two monomials. Writes the result to ec.
function monom_lcm!(ec::ExponentVector{T}, ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    s = zero(T)
    @inbounds for i in 2:length(ec)
        ec[i] = max(ea[i], eb[i])
        s += ec[i]
    end
    ec[1] = s
    monom_overflow_check(ec)
    ec
end

# Checks if the gcd of monomials is constant.
function monom_is_gcd_const(ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    @invariant length(ea) == length(eb)
    @inbounds for j in 2:length(ea)
        if !iszero(ea[j]) && !iszero(eb[j])
            return false
        end
    end
    true
end

# Returns the product of two monomials. Writes the result to ec.
function monom_product!(
    ec::ExponentVector{T},
    ea::ExponentVector{T},
    eb::ExponentVector{T}
) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    @inbounds for j in 1:length(ec)
        ec[j] = ea[j] + eb[j]
    end
    monom_overflow_check(ec)
    ec
end

# Returns the result of monomial division. Writes the result to ec.
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

# Checks monomial divisibility.
function monom_is_divisible(ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    @invariant length(ea) == length(eb)
    @inbounds for j in 2:length(ea)
        if ea[j] < eb[j]
            return false
        end
    end
    return true
end

# Checks monomial divisibility and performs division
function monom_is_divisible!(
    ec::ExponentVector{T},
    ea::ExponentVector{T},
    eb::ExponentVector{T}
) where {T}
    @invariant length(ec) == length(ea) == length(eb)
    !monom_is_divisible(ea, eb) && return false, ec
    monom_division!(ec, ea, eb)
    true, ec
end

# Checks monomials for element-wise equality
function monom_is_equal(ea::ExponentVector{T}, eb::ExponentVector{T}) where {T}
    ea == eb
end

###
# Monomial division masks.

# Constructs the division mask for the given monomial.
#
# If compressed is true, the mask will be compressed so that more variables can
# be handled. Compressed and non-compressed division masks are not compatible.
function monom_create_divmask(
    e::ExponentVector{T},
    _::Type{Mask},
    ndivvars::Int,
    divmap::Vector{U},
    ndivbits::Int,
    compressed::Bool
) where {T, Mask, U}
    ctr = one(Mask)
    res = zero(Mask)
    o = one(Mask)
    if !compressed
        # if the number of variables is small (<= 32)
        @inbounds for i in 1:ndivvars
            for _ in 1:ndivbits
                if e[i + 1] >= divmap[ctr]
                    res |= o << (ctr - 1)
                end
                ctr += o
            end
        end
    else
        # if the number of variables is large (> 32)
        @invariant ndivbits == 1
        @invariant length(divmap) == 8 * sizeof(Mask)
        @inbounds for i in 1:length(divmap)
            nvars_per_bit = divmap[i]
            bit = Mask(0)
            for _ in 1:nvars_per_bit
                bit |= Mask(!iszero(e[ctr + 1]))
                ctr += o
            end
            res |= bit << (i - 1)
        end
        res
    end
    res
end
