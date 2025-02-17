# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Backend for auxiliary function: leading_term

# polynomials => polynomials
function leading_term0(polynomials, options::KeywordArguments)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    result_monoms, result_coeffs = _leading_term1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, result_monoms, result_coeffs, options)
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function leading_term1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _leading_term1(ring, monoms, coeffs, options)
end

function _leading_term1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    params = AlgorithmParameters(ring, options; hint=:large_exponents)
    _, _ring, _monoms, _coeffs = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    _result_monoms, _result_coeffs = leading_term2(_ring, _monoms, _coeffs, params)
    result_monoms, result_coeffs =
        ir_convert_internal_to_ir(_ring, _result_monoms, _result_coeffs, params)
    result_monoms, result_coeffs
end

# internal structs => internal structs
function leading_term2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    result_monoms = map(f -> isempty(f) ? f : [first(f)], monoms)
    result_coeffs = map(f -> isempty(f) ? f : [first(f)], coeffs)
    result_monoms, result_coeffs
end

###
# Backend for auxiliary function: leading_ideal

# polynomials => polynomials
function leading_ideal0(polynomials, options::KeywordArguments)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    result_monoms, result_coeffs = _leading_ideal1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, result_monoms, result_coeffs, options)
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function leading_ideal1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _leading_ideal1(ring, monoms, coeffs, options)
end

function _leading_ideal1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    params = AlgorithmParameters(ring, options; hint=:large_exponents)
    _, _ring, _monoms, _coeffs = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    _result_monoms, _result_coeffs = leading_ideal2(_ring, _monoms, _coeffs, params)
    result_monoms, result_coeffs =
        ir_convert_internal_to_ir(_ring, _result_monoms, _result_coeffs, params)
    result_monoms, result_coeffs
end

# internal structs => internal structs
function leading_ideal2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    if !isgroebner2(ring, monoms, coeffs, params)
        monoms, coeffs = groebner2(ring, monoms, coeffs, params)
    end
    monoms, coeffs = autoreduce2(ring, monoms, coeffs, params)
    result_monoms = map(f -> isempty(f) ? f : [first(f)], monoms)
    result_coeffs = map(f -> isempty(f) ? f : [first(f)], coeffs)
    result_monoms, result_coeffs
end

###
# Backend for auxiliary function: quotient_basis
# Copied from StructuralIdentifiability.jl. The license is MIT.
# Copied from RationalUnivariateRepresentation.jl. The license is MIT.

# polynomials => polynomials
function quotient_basis0(polynomials, options::KeywordArguments)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    result_monoms, result_coeffs = _quotient_basis1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, result_monoms, result_coeffs, options)
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function quotient_basis1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _quotient_basis1(ring, monoms, coeffs, options)
end

function _quotient_basis1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    params = AlgorithmParameters(ring, options; hint=:large_exponents)
    _, _ring, _monoms, _coeffs = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    _result_monoms, _result_coeffs = quotient_basis2(_ring, _monoms, _coeffs, params)
    result_monoms, result_coeffs =
        ir_convert_internal_to_ir(_ring, _result_monoms, _result_coeffs, params)
    result_monoms, result_coeffs
end

# internal structs => internal structs
function quotient_basis2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    monoms, coeffs = leading_ideal2(ring, monoms, coeffs, params)
    # GB = {0}
    if isempty(monoms[1])
        throw(DomainError("Input does not define a zero-dimensional ideal."))
    end
    n = ring.nvars
    leading_exponents = map(m -> monom_to_vector!(zeros(Int, n), first(m)), monoms)
    # GB != {1} and GB is not zero dimensional
    if monoms[1][1] != monom_construct_const(M, n) &&
       length(filter(e -> count(iszero, e) == n - 1, leading_exponents)) < n
        throw(DomainError("Input does not define a zero-dimensional ideal."))
    end
    @invariant dimension2(ring, monoms, coeffs, params) in (-1, 0)
    exponents_to_check = Vector{Vector{Int}}()
    exponents_checked = Set{Vector{Int}}()
    basis_exponents = Set{Vector{Int}}()
    push!(exponents_to_check, zeros(Int, n))
    function _divisible(e1, e2)
        for i in 1:length(e1)
            e1[i] < e2[i] && return false
        end
        true
    end
    while length(exponents_to_check) > 0
        e = popfirst!(exponents_to_check)
        push!(exponents_checked, e)
        if !any(le -> _divisible(e, le), leading_exponents)
            push!(basis_exponents, e)
            for i in 1:n
                next_e = copy(e)
                next_e[i] += 1
                if !(next_e in exponents_checked) && !(next_e in exponents_to_check)
                    push!(exponents_to_check, next_e)
                end
            end
        end
    end
    qb_len = length(basis_exponents)
    basis_exponents = [[monom_construct_from_vector(M, v)] for v in basis_exponents]
    basis_coeffs = Vector{Vector{C}}([[one(coeffs[1][1])] for _ in 1:qb_len])
    return basis_exponents, basis_coeffs
end

###
# Backend for auxiliary function: dimension
# Copied from AlgebraicSolving.jl. The license is GNU GPL, Version 3.0+.

# polynomials => polynomials
function dimension0(polynomials, options::KeywordArguments)
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    _dimension1(ring, monoms, coeffs, options)
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function dimension1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
    _dimension1(ring, monoms, coeffs, options)
end

function _dimension1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    params = AlgorithmParameters(ring, options; hint=:large_exponents)
    _, _ring, _monoms, _coeffs = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    result = dimension2(_ring, _monoms, _coeffs, params)
    result
end

# internal structs => internal structs
function dimension2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    params = deepcopy(params)
    params.target_ord = DegRevLex()
    monoms, coeffs = leading_ideal2(ring, monoms, coeffs, params)
    # GB = {0} or GB = {1}
    if isempty(monoms[1])
        return -1
    end
    res = [trues(ring.nvars)]
    lead_exps = map(x -> monom_to_vector!(zeros(Int, ring.nvars), first(x)), monoms)
    for lexp in lead_exps
        to_del = Int[]
        new_miss = BitVector[]
        for (i, mis) in enumerate(res)
            nz_exps_inds = findall(e -> !iszero(e), lexp)
            ind_var_inds = findall(mis)
            if issubset(nz_exps_inds, ind_var_inds)
                for j in nz_exps_inds
                    new_mis = copy(mis)
                    new_mis[j] = false
                    push!(new_miss, new_mis)
                end
                push!(to_del, i)
            end
        end
        deleteat!(res, to_del)
        append!(res, new_miss)
        unique!(res)
    end
    length(res) == 0 && return -1
    max_length = maximum(mis -> length(findall(mis)), res)
    return max_length
end
