# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Input-output conversion of polynomials.

# Our conventions on some edge cases:
# - Trying to compute a Groebner basis of an empty set is an error.
# - The Groebner basis of [0] is [0].
# - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]

# Polynomials from frontend packages, such as AbstractAlgebra.jl and
# DynamicPolynomials.jl, may not support some of the orderings supported by
# Groebner.jl. We compute the basis in the requested ordering in Groebner.jl,
# and then order the terms in the output according to the ordering of the input
# polynomials.

# Zero polynomial is represented with an empty vector of monomials and an empty
# vector of coefficients.

@noinline function __throw_input_not_supported(msg, val)
    throw(DomainError(val, "This input is not supported, sorry.\n$msg"))
end

###
# Converting frontend polynomials to intermediate representation (ir)

# Polynomials in the intermediate representation are represented 
# by coefficients (UInt64) and exponent vectors (Vector{UInt64}).

"""
    PolyRing

Polynomial ring.

## Example

4 variables, modulo 2^30 + 3, degrevlex order.

```julia
ring = Groebner.PolyRing(4, Groebner.DegRevLex(), 2^30+3)
```

4 variables, the rationals, a block order.

```julia
ord = Groebner.DegRevLex([1,2]) * Groebner.DegRevLex([3,4])
ring = Groebner.PolyRing(4, ord, 0)
```
"""
mutable struct PolyRing{
    Ord <: AbstractMonomialOrdering,
    C <: Union{CoeffZp, CompositeCoeffZp}
}
    nvars::Int
    ord::Ord
    ch::C
end

io_iszero_coeffs(v) = isempty(v)
io_iszero_monoms(v) = isempty(v)

io_zero_coeffs(::Type{T}, ring::PolyRing) where {T} = Vector{T}()
io_zero_monoms(::Type{T}, ring::PolyRing) where {T} = Vector{T}()

function ir_basic_is_valid(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}}
) where {T <: Integer, C <: Number}
    !(length(monoms) == length(coeffs)) && throw(DomainError("Invalid IR."))
    isempty(monoms) && throw(DomainError("Invalid IR."))
    !(ring.nvars > 0) && throw(DomainError("The number of variables must be positive."))
    !(ring.ch >= 0) && throw(DomainError("Field characteristic must be nonnegative."))
    if ring.ch > 0
        !(C <: Integer) && throw(DomainError("Coefficients must be integers."))
        !(ring.ch <= typemax(C)) && throw(DomainError("Invalid IR."))
    else
        !(C <: Rational || C <: Integer) &&
            throw(DomainError("Coefficients must be integer or rationals."))
    end
    (ring.ord == InputOrdering()) && throw(DomainError("Invalid IR."))
end

function ir_is_valid(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}}
) where {T <: Integer, C <: Number}
    ir_basic_is_valid(ring, monoms, coeffs)
    for i in 1:length(monoms)
        !(length(monoms[i]) == length(coeffs[i])) && throw(DomainError("Invalid IR."))
        for j in 1:length(monoms[i])
            !(length(monoms[i][j]) == ring.nvars) && throw(DomainError("Invalid IR."))
            !(all(>=(0), monoms[i][j])) && throw(DomainError("Invalid IR."))
            iszero(coeffs[i][j]) && throw(DomainError("Invalid IR")) # can be relaxed
            if (ring.ch > 0)
                !(0 < coeffs[i][j] < ring.ch) && throw(DomainError("Invalid IR."))
            end
        end
    end
    true
end

ir_ensure_assumptions(ring, monoms, coeffs) =
    throw(DomainError("Invalid IR, unknown types."))

function ir_ensure_assumptions(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{M}}},
    coeffs::Vector{Vector{C}}
) where {M <: Integer, C <: Coeff}
    ir_basic_is_valid(ring, monoms, coeffs)
    # Copy input
    new_monoms, new_coeffs = empty(monoms), empty(coeffs)
    for i in 1:length(monoms)
        !(length(monoms[i]) == length(coeffs[i])) && throw(DomainError("Invalid IR."))
        push!(new_monoms, empty(monoms[i]))
        push!(new_coeffs, empty(coeffs[i]))
        for j in 1:length(monoms[i])
            !(length(monoms[i][j]) == ring.nvars) && throw(DomainError("Invalid IR."))
            !(all(>=(0), monoms[i][j])) && throw(DomainError("Invalid IR."))
            push!(new_monoms[i], monoms[i][j])
            push!(new_coeffs[i], coeffs[i][j])
        end
    end
    # Normalize
    for i in 1:length(new_monoms)
        for j in 1:length(new_monoms[i])
            if ring.ch > 0
                new_coeffs[i][j] = mod(new_coeffs[i][j], ring.ch)
            end
        end
    end
    # Sort
    # if !issorted(new_monoms, by=(i, j) -> monom_isless(i, j, ring.ord))
    #     perm = collect(1:length(new_monoms))
    #     sort!(perm, by=(i, j) -> monom_isless(new_monoms[j], new_monoms[i], ring.ord))
    #     new_monoms = new_monoms[perm]
    #     new_coeffs = new_coeffs[perm]
    # end
    # Normalize
    _new_monoms = empty(new_monoms)
    _new_coeffs = empty(new_coeffs)
    for i in 1:length(new_monoms)
        push!(_new_monoms, empty(new_monoms[i]))
        push!(_new_coeffs, empty(new_coeffs[i]))
        if !isempty(new_coeffs[i])
            push!(_new_monoms[i], new_monoms[i][1])
            push!(_new_coeffs[i], new_coeffs[i][1])
        end
        slow_idx = 1
        for fast_idx in 2:length(new_monoms[i])
            if _new_monoms[i][slow_idx] != new_monoms[i][fast_idx]
                push!(_new_monoms[i], new_monoms[i][fast_idx])
                push!(_new_coeffs[i], new_coeffs[i][fast_idx])
                slow_idx += 1
                continue
            end
            _new_coeffs[i][slow_idx] =
                Base.Checked.checked_add(_new_coeffs[i][slow_idx], new_coeffs[i][fast_idx])
            if ring.ch > 0 && _new_coeffs[i][slow_idx] >= ring.ch
                _new_coeffs[i][slow_idx] -= ring.ch
                @invariant _new_coeffs[i][slow_idx] < ring.ch
            end
        end
    end
    new_monoms, new_coeffs = _new_monoms, _new_coeffs
    # Remove zeros
    for i in 1:length(new_monoms)
        perm = collect(1:length(new_monoms[i]))
        filter!(j -> !iszero(new_coeffs[i][j]), perm)
        new_monoms[i] = new_monoms[i][perm]
        new_coeffs[i] = new_coeffs[i][perm]
    end
    ring, new_monoms, new_coeffs
end

function io_convert_polynomials_to_ir(polynomials, options::KeywordArguments)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring = io_extract_ring(polynomials)
    var_to_index, monoms, coeffs = _io_convert_polynomials_to_ir(ring, polynomials)
    ring = PolyRing(ring.nvars, ordering_transform(ring.ord, var_to_index), ring.ch)
    options.ordering = ordering_transform(options.ordering, var_to_index)
    ring, monoms, coeffs, options
end

function _io_convert_polynomials_to_ir(ring::PolyRing, polynomials)
    coeffs = io_extract_coeffs_ir(ring, polynomials)
    reversed_order, var_to_index, monoms = io_extract_monoms_ir(ring, polynomials)
    @invariant length(coeffs) == length(monoms)
    var_to_index, monoms, coeffs
end

function io_convert_ir_to_polynomials(
    ring::PolyRing,
    polynomials,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    options
) where {M <: Monom, C <: Coeff}
    @assert !isempty(polynomials)
    _io_convert_ir_to_polynomials(ring, polynomials, monoms, coeffs, options)
end

###
# Converting to internal representation

function io_convert_ir_to_internal(ring, monoms, coeffs, params, repr)
    monoms2 = Vector{Vector{repr.monomtype}}(undef, length(monoms))
    coeffs2 = Vector{Vector{repr.coefftype}}(undef, length(monoms))
    @inbounds for i in 1:length(monoms)
        monoms2[i] = Vector{repr.monomtype}(undef, length(monoms[i]))
        coeffs2[i] = Vector{repr.coefftype}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            monoms2[i][j] = monom_construct_from_vector(repr.monomtype, monoms[i][j])
            coeffs2[i][j] = repr.coefftype(coeffs[i][j])
        end
    end
    ring2, term_sorting_permutations =
        io_set_monomial_ordering!(ring, monoms2, coeffs2, params)
    term_sorting_permutations, ring2, monoms2, coeffs2
end

function io_convert_internal_to_ir(ring, monoms, coeffs, params)
    monoms2 = Vector{Vector{Vector{UInt64}}}(undef, length(monoms))
    coeffs2 = coeffs
    @inbounds for i in 1:length(monoms)
        monoms2[i] = Vector{Vector{UInt64}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            monoms2[i][j] = Vector{UInt64}(undef, ring.nvars)
            monom_to_vector!(monoms2[i][j], monoms[i][j])
        end
    end
    monoms2, coeffs2
end

###
# Converting to selected polynomial representation

function io_convert_to_internal(
    representation::PolynomialRepresentation,
    polynomials,
    kws::KeywordArguments;
    dropzeros=true
)
    io_check_input(polynomials, kws)
    # NOTE: Input polynomials must not be modified.
    @log :misc "Converting input polynomials to internal representation.."
    ring = io_extract_ring(polynomials)
    var_to_index, monoms, coeffs =
        io_extract_polynomial_data(representation, ring, polynomials)
    @log :misc "Done converting input polynomials to internal representation."
    @log :all """
    Polynomials in internal representation:
    Ring: $ring
    Variable to index map: $var_to_index
    Monomials: $monoms
    Coefficients: $coeffs"""
    if dropzeros
        @log :misc "Removing zero polynomials"
        io_remove_zeros_from_input!(ring, monoms, coeffs)
    end
    ring, var_to_index, monoms, coeffs
end

function io_check_input(polynomials, kws)
    if isempty(polynomials)
        __throw_input_not_supported(polynomials, "Empty input array.")
    end
    _io_check_input(polynomials, kws)
end

function io_extract_polynomial_data(
    representation::PolynomialRepresentation,
    ring::PolyRing,
    polynomials::Vector{T}
) where {T}
    coeffs = io_extract_coeffs(representation, ring, polynomials)
    reversed_order, var_to_index, monoms =
        io_extract_monoms(representation, ring, polynomials)
    @invariant length(coeffs) == length(monoms)
    if reversed_order
        for i in 1:length(coeffs)
            @invariant length(coeffs[i]) == length(monoms[i])
            reverse!(coeffs[i])
            reverse!(monoms[i])
        end
    end
    var_to_index, monoms, coeffs
end

function io_remove_zeros_from_input!(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{T}}
) where {M, T}
    @invariant length(monoms) == length(coeffs)
    filter!(!io_iszero_coeffs, coeffs)
    filter!(!io_iszero_monoms, monoms)
    @invariant length(monoms) == length(coeffs)
    iszerobasis = isempty(monoms)
    @log :all "After removing zero polynomials from input:" monoms coeffs
    iszerobasis
end

# Checks that the monomial orderings specified by the given `ring` and
# `params.target_ord` are consistent with the given input monomials `monoms`. In
# case the target ordering differs from the `ring` ordering,  
# sorts the polynomials terms w.r.t. the target ordering.
#
# Also returns the sorting permutations for polynomial terms
function io_set_monomial_ordering!(ring, monoms, coeffs, params)
    if ring.ord == params.target_ord
        # No reordering of terms needed, the terms are already ordered according
        # to the requested monomial ordering
        return ring, Vector{Vector{Int}}()
    end
    ring = PolyRing(ring.nvars, params.target_ord, ring.ch)
    permutations = sort_input_terms_to_change_ordering!(monoms, coeffs, params.target_ord)
    ring, permutations
end

###
# Converting polynomials from internal representation back to original types

function io_convert_to_output(
    ring::PolyRing,
    polynomials,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @assert !isempty(polynomials)
    @log :misc "Converting polynomials from internal representation to output format"
    # NOTE: Internal polynomials must not be modified.
    if isempty(monoms)
        @log :debug "Output is empty, appending an empty placeholder polynomial"
        push!(monoms, io_zero_monoms(M, ring))
        push!(coeffs, io_zero_coeffs(C, ring))
    end
    _io_convert_to_output(ring, polynomials, monoms, coeffs, params)
end

function io_convert_changematrix_to_output(
    ring::PolyRing,
    polynomials,
    npolys::Int,
    monoms::Vector{Vector{Vector{M}}},
    coeffs::Vector{Vector{Vector{C}}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    @assert !isempty(polynomials)
    @log :misc "Converting polynomials from internal representation to output format"
    changematrix = Matrix{eltype(polynomials)}(undef, length(monoms), length(polynomials))
    for i in 1:length(monoms)
        matrix_row = io_convert_to_output(ring, polynomials, monoms[i], coeffs[i], params)
        matrix_row_full = Vector{eltype(matrix_row)}(undef, length(polynomials))
        k = 1
        for j in 1:length(polynomials)
            if iszero(polynomials[j])
                matrix_row_full[j] = zero(polynomials[j])
            else
                matrix_row_full[j] = matrix_row[k]
                k += 1
            end
        end
        changematrix[i, :] .= matrix_row_full
    end
    changematrix
end

###
# Utilities for composite coefficients

function io_convert_to_output_batched(
    ring::PolyRing,
    batch::NTuple{N, V},
    monoms::Vector{Vector{M}},
    coeffs_packed::Vector{Vector{C}},
    params::AlgorithmParameters
) where {N, V, M <: Monom, C <: Coeff}
    @assert !any(isempty, batch)
    @log :misc "Converting polynomials from internal representation to output format (batched)"

    coeffs_unpacked = io_unpack_composite_coefficients(coeffs_packed)
    # Internal polynomials must not be modified.
    if isempty(monoms)
        @log :debug "Output is empty, appending an empty placeholder polynomial"
        push!(monoms, io_zero_monoms(M, ring))
        for coeffs in coeffs_unpacked
            push!(coeffs, io_zero_coeffs(C, ring))
        end
    end

    ntuple(i -> io_convert_to_output(ring, batch[i], monoms, coeffs_unpacked[i], params), N)
end

function io_unpack_composite_coefficients(
    composite_coeffs::Vector{Vector{CompositeNumber{N, T}}}
) where {N, T <: CoeffZp}
    coeffs_part_i = ntuple(_ -> Vector{Vector{T}}(undef, length(composite_coeffs)), N)
    @inbounds for i in 1:length(composite_coeffs)
        for k in 1:N
            coeffs_part_i[k][i] = Vector{T}(undef, length(composite_coeffs[i]))
        end
        for j in 1:length(composite_coeffs[i])
            ai = composite_coeffs[i][j].data
            for k in 1:N
                coeffs_part_i[k][i][j] = ai[k]
            end
        end
    end
    coeffs_part_i
end
