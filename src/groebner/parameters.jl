# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Select parameters in Groebner basis computation

# Specifies linear algebra algorithm
struct LinearAlgebra
    algorithm::Symbol
    sparsity::Symbol
end

# Specifies polynomial representation
struct PolynomialRepresentation
    monomtype::Any
    coefftype::Any
    # If this field is false, then any implementation of the arithmetic in Z/Zp
    # must cast the coefficients into a wider integer type before performing any
    # arithmetic operations to avoid the risk of overflow.
    using_wide_type_for_coeffs::Bool
end

function param_select_polynomial_representation(
    char::Coeff,
    nvars::Int,
    ground::Symbol,
    ordering::AbstractMonomialOrdering,
    homogenize::Bool,
    monoms::Symbol,
    arithmetic::Symbol;
    hint::Symbol=:none
)
    if !(hint in (:none, :large_exponents))
        @info "The given hint=$hint was discarded"
    end
    monomtype = param_select_monomtype(char, nvars, ordering, homogenize, hint, monoms)
    coefftype, using_wide_type_for_coeffs =
        param_select_coefftype(char, nvars, ground, ordering, homogenize, hint, monoms, arithmetic)
    PolynomialRepresentation(monomtype, coefftype, using_wide_type_for_coeffs)
end

function param_select_monomtype(
    char::Coeff,
    nvars::Int,
    ordering::AbstractMonomialOrdering,
    homogenize::Bool,
    hint::Symbol,
    monoms::Symbol
)
    if hint === :large_exponents
        # use 64 bits if large exponents detected
        desired_monom_type = ExponentVector{UInt64}
        @assert monom_is_supported_ordering(desired_monom_type, ordering)
        return desired_monom_type
    end

    # If homogenization is requested, or if a part of the ordering is
    # lexicographical, the generators will potentially be homogenized later.
    if homogenize
        desired_monom_type = ExponentVector{UInt32}
        @assert monom_is_supported_ordering(desired_monom_type, ordering)
        return desired_monom_type
    end

    ExponentSize = UInt8
    variables_per_word = div(sizeof(UInt), sizeof(ExponentSize))
    # if dense representation is requested
    if monoms === :dense
        @assert monom_is_supported_ordering(ExponentVector{ExponentSize}, ordering)
        return ExponentVector{ExponentSize}
    end

    # if packed representation is requested
    if monoms === :packed
        if monom_is_supported_ordering(PackedTuple1{UInt64, ExponentSize}, ordering)
            if nvars < variables_per_word
                return PackedTuple1{UInt64, ExponentSize}
            elseif nvars < 2 * variables_per_word
                return PackedTuple2{UInt64, ExponentSize}
            elseif nvars < 3 * variables_per_word
                return PackedTuple3{UInt64, ExponentSize}
            elseif nvars < 4 * variables_per_word
                return PackedTuple4{UInt64, ExponentSize}
            end
            # falling back to dense representation
        end
    end

    # in the automatic choice, we always prefer packed representations
    if monoms === :auto
        if monom_is_supported_ordering(PackedTuple1{UInt64, ExponentSize}, ordering)
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

tight_signed_int_type(x::CompositeNumber{N, T}) where {N, T} =
    CompositeNumber{N, mapreduce(tight_signed_int_type, promote_type, x.data)}
tight_unsigned_int_type(x::CompositeNumber{N, T}) where {N, T} =
    CompositeNumber{N, mapreduce(tight_signed_int_type, promote_type, x.data)}

function tight_signed_int_type(x::T) where {T <: Integer}
    types = (Int8, Int16, Int32, Int64, Int128)
    idx = findfirst(T -> x <= typemax(T), types)
    @assert !isnothing(idx)
    types[idx]
end

function tight_unsigned_int_type(x::T) where {T <: Integer}
    types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    idx = findfirst(T -> x <= typemax(T), types)
    @assert !isnothing(idx)
    types[idx]
end

function param_select_coefftype(
    char::Coeff,
    nvars::Int,
    ground::Symbol,
    ordering::AbstractMonomialOrdering,
    homogenize::Bool,
    hint::Symbol,
    monoms::Symbol,
    arithmetic::Symbol;
    using_wide_type_for_coeffs::Bool=false
)
    if ground == :generic
        return CoeffGeneric, true
    end
    if iszero(char)
        return Rational{BigInt}, true
    end
    @assert char > 0
    @assert char < typemax(UInt64)

    tight_signed_type = tight_signed_int_type(char)

    if arithmetic === :signed
        if typemax(Int32) < char < typemax(UInt32) || typemax(Int64) < char < typemax(UInt64)
            @info "Cannot use $(arithmetic) arithmetic with characteristic $char"
            @assert false
        elseif !using_wide_type_for_coeffs
            return tight_signed_type, using_wide_type_for_coeffs
        else
            return widen(tight_signed_type), using_wide_type_for_coeffs
        end
    end

    tight_unsigned_type = tight_unsigned_int_type(char)
    tight_unsigned_type = if !using_wide_type_for_coeffs
        tight_unsigned_type
    else
        widen(tight_unsigned_type)
    end

    tight_unsigned_type, using_wide_type_for_coeffs
end

# Stores parameters for a single GB computation.
struct AlgorithmParameters{Arithmetic <: AbstractArithmetic}
    # Desired monomial ordering of output polynomials
    target_ord::AbstractMonomialOrdering
    # Original monomial ordering of input polynomials
    original_ord::AbstractMonomialOrdering

    # Specifies correctness checks levels
    heuristic_check::Bool
    randomized_check::Bool
    certify_check::Bool

    # If do homogenize input generators
    homogenize::Bool

    # This option only makes sense for some functions (e.g., `normalform`). It
    # specifies if we should check if the input is indeed a Groebner basis.
    check::Bool

    # Linear algebra backend
    linalg::LinearAlgebra

    # Arithmetic in the ground field
    arithmetic::Arithmetic

    # Representation of polynomials in F4
    representation::PolynomialRepresentation

    # If reduced Groebner basis is needed
    reduced::Bool

    # Strategy for modular computation in groebner. This can be one of the
    # following:
    # - :classic_modular
    # - :learn_and_apply
    modular_strategy::Symbol

    # The width of composite numbers in learn & apply
    composite::Int

    # Multi-threading
    threaded_f4::Symbol
    threaded_multimodular::Symbol

    rng::Random.Xoshiro

    changematrix::Bool
end

function AlgorithmParameters(ring::PolyRing, kwargs::KeywordArguments; hint=:none)
    if kwargs.ordering === InputOrdering() || kwargs.ordering === nothing
        ordering = ring.ord
    else
        ordering = kwargs.ordering
    end
    target_ord = ordering
    original_ord = ring.ord

    heuristic_check = true
    randomized_check = true
    certify_check = kwargs.certify

    homogenize = if kwargs.homogenize === :yes
        true
    else
        if kwargs.homogenize === :auto
            if ring.nvars <= 1
                false
            elseif target_ord isa Lex || target_ord isa ProductOrdering
                true
            else
                false
            end
        else
            false
        end
    end

    # Over Z_p, linalg = :randomized signals randomized linear algebra.
    # Over Q, linalg = :randomized has almost no effect since multi-modular
    # tracing is used by default. Still, we say "almost" since some routines in
    # multi-modular tracing may compute Groebner bases in Zp for
    # checking/verification, and they benefit from linalg = :randomized.
    linalg = kwargs.linalg
    if ring.ground === :zp && (linalg === :randomized || linalg === :auto)
        # Do not use randomized linear algebra if the field characteristic is
        # too small. 
        # TODO: In the future, it would be good to adapt randomized linear
        # algebra to this case by taking more random samples
        if ring.characteristic < 500
            if linalg === :randomized
                @info """
                The option linalg=:randomized was ignored:
                the ground field characteristic is small."""
            end
            linalg = :deterministic
        end
    end
    if ring.ground == :generic
        linalg = :deterministic
    end
    if linalg === :auto
        linalg = :randomized
    end
    if kwargs.function_id === :isgroebner
        linalg =
            if ring.ground === :zp && (kwargs.linalg === :randomized || kwargs.linalg === :auto)
                :randomized
            else
                :deterministic
            end
    end
    linalg_sparsity = :sparse
    linalg_algorithm = LinearAlgebra(linalg, linalg_sparsity)

    representation = param_select_polynomial_representation(
        ring.characteristic,
        ring.nvars,
        ring.ground,
        target_ord,
        homogenize,
        kwargs.monoms,
        kwargs.arithmetic,
        hint=hint
    )

    arithmetic = select_arithmetic(
        representation.coefftype,
        ring.characteristic,
        kwargs.arithmetic,
        representation.using_wide_type_for_coeffs
    )

    reduced = kwargs.reduced

    threaded = kwargs.threaded
    if !(_threaded[])
        if threaded === :yes
            @info """
            Keyword argument `threaded = :yes` was provided to Groebner.jl,
            however, multi-threading is disabled globally in Groebner.jl due to
            the environment variable GROEBNER_NO_THREADED = 0.

            Consider enabling threading by setting GROEBNER_NO_THREADED to 1."""
        end
        threaded = :no
    end

    if ring.ground === :zp
        threaded_f4 = threaded
        threaded_multimodular = :no
    elseif ring.ground == :qq
        threaded_f4 = :no
        threaded_multimodular = threaded
    else
        threaded_f4 = :no
        threaded_multimodular = :no
    end

    # By default, modular computation uses learn & apply
    modular_strategy = kwargs.modular
    if modular_strategy === :auto
        modular_strategy = :learn_and_apply
    end
    if !reduced
        # The option reduced=false was passed in the input, 
        # falling back to classic multi-modular algorithm.
        modular_strategy = :classic_modular
    end
    composite = kwargs._composite

    seed = kwargs.seed
    rng = Random.Xoshiro(seed)

    changematrix = kwargs.changematrix
    if changematrix
        if !(target_ord isa DegRevLex)
            throw(DomainError("Only DegRevLex is supported with changematrix = true."))
        end
        if (ring.ground == :generic)
            throw(DomainError("Generic fields are not supported with changematrix = true."))
        end
    end

    if kwargs.function_id == :groebner_learn || kwargs.function_id == :groebner_apply!
        if ring.ground == :generic
            throw(DomainError("Generic fields are not supported with learn & apply."))
        end
    end

    AlgorithmParameters(
        target_ord,
        original_ord,
        heuristic_check,
        randomized_check,
        certify_check,
        homogenize,
        kwargs.check,
        linalg_algorithm,
        arithmetic,
        representation,
        reduced,
        modular_strategy,
        composite,
        threaded_f4,
        threaded_multimodular,
        rng,
        changematrix
    )
end

function param_mod_p(
    params::AlgorithmParameters,
    prime::C;
    using_wide_type_for_coeffs=nothing
) where {C <: Coeff}
    is_wide_type_coeffs = if !isnothing(using_wide_type_for_coeffs)
        using_wide_type_for_coeffs
    else
        params.representation.using_wide_type_for_coeffs
    end
    representation = PolynomialRepresentation(
        params.representation.monomtype,
        params.representation.coefftype,
        is_wide_type_coeffs
    )
    struct_update(
        AlgorithmParameters,
        params,
        (
            arithmetic=select_arithmetic(C, prime, :auto, is_wide_type_coeffs),
            representation=representation
        )
    )
end
