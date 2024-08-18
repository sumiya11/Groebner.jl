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
    Ord <: Union{AbstractMonomialOrdering, AbstractInternalOrdering},
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

function io_convert_polynomials_to_ir(polynomials, options::KeywordArguments)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring = io_extract_ring(polynomials)
    var_to_index, monoms, coeffs = _io_convert_polynomials_to_ir(ring, polynomials)
    ring = PolyRing(ring.nvars, io_convert_ordering_to_ir(var_to_index, ring.ord), ring.ch)
    options.ordering = io_convert_ordering_to_ir(var_to_index, options.ordering)
    ring, monoms, coeffs
end

function io_convert_ordering_to_ir(var_to_index, ordering)
    Ords = (Lex, DegLex, DegRevLex)
    idx = findfirst(Ord -> ordering isa Ord, Ords)
    if !isnothing(idx)
        if isnothing(ordering.variables)
            Ords[idx]()
        else
            Ords[idx](map(k -> var_to_index[k], ordering.variables))
        end
    elseif ordering isa ProductOrdering
        ProductOrdering(
            io_convert_ordering_to_ir(var_to_index, ordering.ord1),
            io_convert_ordering_to_ir(var_to_index, ordering.ord2)
        )
    else
        ordering
    end
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
    ring2, term_sorting_permutations = io_set_monomial_ordering!(
        ring,
        Dict{Int, Int}(collect(1:(ring.nvars)) .=> collect(1:(ring.nvars))),
        monoms2,
        coeffs2,
        params
    )
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
# Selecting polynomial representation

struct PolynomialRepresentation
    monomtype::Type
    coefftype::Type
    # If this field is false, then any implementation of the arithmetic in Z/Zp
    # must cast the coefficients into a wider integer type before performing any
    # arithmetic operations to avoid the risk of overflow.
    using_wide_type_for_coeffs::Bool
end

io_peek_at_polynomials(polynomials::AbstractVector) =
    __throw_input_not_supported("", polynomials)

function io_select_polynomial_representation(
    ring::PolyRing,
    kws::KeywordArguments;
    hint::Symbol=:none
)
    # isempty(polynomials) && __throw_input_not_supported("Empty input.", polynomials)
    if !(hint in (:none, :large_exponents))
        @log :warn "The given hint=$hint was discarded"
    end
    # frontend, npolys, char, nvars, ordering = io_peek_at_polynomials(polynomials)
    frontend, npolys, char, nvars, ordering =
        :abstractalgebra, 1, ring.ch, ring.nvars, ring.ord
    monomtype = io_select_monomtype(char, nvars, ordering, kws, hint)
    coefftype, using_wide_type_for_coeffs =
        io_select_coefftype(char, nvars, ordering, kws, hint)
    basering = iszero(char) ? :qq : :zp
    @log :misc "Frontend: $frontend"
    @log :misc """
    Input: $npolys polynomials over $basering in $nvars variables
    Ordering: $ordering"""
    @log :misc """
    Internal representation: 
    monomials are $monomtype
    coefficients are $coefftype
    wide type for coefficients: $using_wide_type_for_coeffs"""
    PolynomialRepresentation(monomtype, coefftype, using_wide_type_for_coeffs)
end

function io_select_monomtype(char, nvars, ordering, kws, hint)
    @log :misc """
    Selecting monomial representation.
    Given hint hint=$hint. Keyword argument monoms=$(kws.monoms)"""
    if hint === :large_exponents
        @log :misc "As hint=$hint was provided, using 64 bits per single exponent"
        # use 64 bits if large exponents detected
        desired_monom_type = ExponentVector{UInt64}
        @assert monom_is_supported_ordering(desired_monom_type, kws.ordering)
        return desired_monom_type
    end

    # If homogenization is requested, or if a part of the ordering is
    # lexicographical, the ideal will potentially be homogenized later.
    if kws.homogenize === :yes || (
        kws.homogenize === :auto && (
            kws.ordering isa Lex ||
            kws.ordering isa ProductOrdering ||
            (
                (ordering isa Lex || ordering isa ProductOrdering) &&
                kws.ordering isa InputOrdering
            )
        )
    )
        @log :misc """
        As homogenize=:yes/:auto was provided, 
        representing monomials as dense exponent vectors"""
        desired_monom_type = ExponentVector{UInt32}
        @assert monom_is_supported_ordering(desired_monom_type, kws.ordering)
        return desired_monom_type
    end

    ExponentSize = UInt8
    variables_per_word = div(sizeof(UInt), sizeof(ExponentSize))
    # if dense representation is requested
    if kws.monoms === :dense
        @assert monom_is_supported_ordering(ExponentVector{ExponentSize}, kws.ordering)
        return ExponentVector{ExponentSize}
    end

    # if sparse representation is requested
    if kws.monoms === :sparse
        if monom_is_supported_ordering(
            SparseExponentVector{ExponentSize, Int32, nvars},
            kws.ordering
        )
            return SparseExponentVector{ExponentSize, Int32, nvars}
        end
        @log :info """
        The given monomial ordering $(kws.ordering) is not implemented for
        $(kws.monoms) monomial representation. Falling back to other monomial
        representations."""
    end

    # if packed representation is requested
    if kws.monoms === :packed
        if monom_is_supported_ordering(PackedTuple1{UInt64, ExponentSize}, kws.ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, ExponentSize}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, ExponentSize}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, ExponentSize}
            elseif nvars < 4 * variables_per_word
                return PackedTuple4{UInt64, ExponentSize}
            end
            @log :info """
            Unable to use $(kws.monoms) monomial representation, too many
            variables ($nvars). Falling back to dense monomial
            representation."""
        else
            @log :info """
            The given monomial ordering $(kws.ordering) is not implemented for
            $(kws.monoms) monomial representation. Falling back to dense
            representation."""
        end
    end

    # in the automatic choice, we always prefer packed representations
    if kws.monoms === :auto
        # TODO: also check that ring.ord is supported
        if monom_is_supported_ordering(PackedTuple1{UInt64, ExponentSize}, kws.ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, ExponentSize}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, ExponentSize}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, ExponentSize}
            elseif nvars < 4 * variables_per_word
                return PackedTuple4{UInt64, ExponentSize}
            end
        end
    end

    ExponentVector{ExponentSize}
end

function io_get_tight_signed_int_type(x::T) where {T <: Integer}
    types = (Int8, Int16, Int32, Int64, Int128)
    idx = findfirst(T -> x <= typemax(T), types)
    @assert !isnothing(idx)
    types[idx]
end

function io_get_tight_unsigned_int_type(x::T) where {T <: Integer}
    types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    idx = findfirst(T -> x <= typemax(T), types)
    @assert !isnothing(idx)
    types[idx]
end

function io_select_coefftype(
    char,
    nvars,
    ordering,
    kws,
    hint;
    using_wide_type_for_coeffs=false
)
    @log :misc """
    Selecting coefficient representation.
    Given hint hint=$hint. Keyword argument arithmetic=$(kws.arithmetic)"""

    if iszero(char)
        return Rational{BigInt}, true
    end
    @assert char > 0

    if char > typemax(UInt64)
        __throw_input_not_supported(
            char,
            "The coefficient field characteristic is too large."
        )
    end

    tight_signed_type = io_get_tight_signed_int_type(char)

    if kws.arithmetic === :floating
        return Float64, true
    end

    if kws.arithmetic === :signed
        if typemax(Int32) < char < typemax(UInt32) ||
           typemax(Int64) < char < typemax(UInt64)
            @log :warn "Cannot use $(kws.arithmetic) arithmetic with characteristic $char"
            @assert false
        elseif !using_wide_type_for_coeffs
            return tight_signed_type, using_wide_type_for_coeffs
        else
            return widen(tight_signed_type), using_wide_type_for_coeffs
        end
    end

    tight_unsigned_type = io_get_tight_unsigned_int_type(char)
    tight_unsigned_type = if !using_wide_type_for_coeffs
        tight_unsigned_type
    else
        widen(tight_unsigned_type)
    end

    tight_unsigned_type, using_wide_type_for_coeffs
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
function io_set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)
    current_ord = ring.ord
    target_ord = params.target_ord
    internal_ord = io_convert_to_internal_monomial_ordering(var_to_index, target_ord)
    @log :misc "Internal ordering:\n$internal_ord"
    ring = PolyRing(ring.nvars, internal_ord, ring.ch)
    if current_ord == target_ord
        # No reordering of terms needed, the terms are already ordered according
        # to the requested monomial ordering
        return ring, Vector{Vector{Int}}()
    end
    @log :misc "Reordering input polynomial terms from $(current_ord) to $(target_ord)"
    permutations = sort_input_terms_to_change_ordering!(monoms, coeffs, internal_ord)
    @log :all "Reordered terms:" monoms coeffs
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
