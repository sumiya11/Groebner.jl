# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Intermediate representation of polynomials (ir)

# Polynomials in the intermediate representation are represented by two arrays:
# an array of exponent vectors and an array of coefficients. Zero polynomial is
# represented with two empty arrays.

const IR_exponent = UInt32
const IR_coeff = Union{Number, CoeffGeneric}

# Intermediate representation (ir) enforces a number assumtpions:
# - Intermediate representation owns the memory.
# - Monomials in intermediate representation are sorted and unique.
# - Coefficients in intermediate representation are nonzero and normalized.

# Our conventions on some edge cases:
# - Trying to compute a Groebner basis of an empty set is an error.
# - The Groebner basis of [0] is [0].
# - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]

struct PolyRing{Ord <: AbstractMonomialOrdering, C <: Union{CoeffZp, CompositeCoeffZp}}
    nvars::Int
    ord::Ord
    characteristic::C
    ground::Symbol

    PolyRing(nvars::Int, ord, ch) = PolyRing(nvars, ord, ch, iszero(ch) ? :qq : :zp)

    function PolyRing(nvars::Int, ord::Ord, ch::C, ground::Symbol) where {Ord, C}
        !(ground in (:zp, :qq, :generic)) && throw(DomainError("Invalid ground field."))
        !(nvars >= 0) && throw(DomainError("The number of variables must be non-negative."))
        new{Ord, C}(nvars, ord, ch, ground)
    end
end

Base.:(==)(r1::PolyRing, r2::PolyRing) =
    r1.nvars == r2.nvars &&
    r1.ord == r2.ord &&
    r1.characteristic == r2.characteristic &&
    r1.ground == r2.ground

ir_is_valid_basic(batch) = throw(DomainError("Invalid IR, unknown types."))
ir_is_valid_basic(ring, monoms, coeffs) = throw(DomainError("Invalid IR, unknown types."))

function ir_is_valid_basic(batch::NTuple{N, T}) where {N, T}
    for el in batch
        length(el) != 3 && throw(DomainError("Invalid IR."))
    end
    for el in batch
        ir_is_valid_basic(el...)
    end
end

function ir_is_valid_basic(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{E}}},
    coeffs::Vector{Vector{C}}
) where {E <: Integer, C <: IR_coeff}
    !(length(monoms) == length(coeffs)) && throw(DomainError("Invalid IR."))
    isempty(monoms) && throw(DomainError("Invalid IR."))
    !(ring.nvars >= 0) && throw(DomainError("The number of variables must be non-negative."))
    !(ring.characteristic >= 0) && throw(DomainError("Field characteristic must be nonnegative."))
    if ring.ground == :zp
        !(C <: Integer) && throw(DomainError("Coefficients must be integers."))
        (C <: BigInt) && throw(DomainError("Coefficients must fit in a machine register."))
        !(ring.characteristic <= typemax(C)) && throw(DomainError("Invalid IR."))
        !Primes.isprime(ring.characteristic) &&
            throw(DomainError("Ring characteristic must be prime"))
    elseif ring.ground == :qq
        !(C <: Rational || C <: Integer || C <: CoeffGeneric) &&
            throw(DomainError("Coefficients must be integer or rational."))
    else
        !(C <: CoeffGeneric) && throw(DomainError("Coefficients must be CoeffGeneric."))
    end
    (ring.ord == InputOrdering()) && throw(DomainError("Invalid monomial ordering."))
    vars = ordering_variables(ring.ord)
    if !isempty(vars) && !(Set(vars) == Set(1:(ring.nvars)))
        throw(
            DomainError(
                "Invalid monomial ordering. Expected variables to be numbers from 1 to $(ring.nvars), but got: $(vars)"
            )
        )
    end
end

function ir_is_valid(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}}
) where {I <: Integer, C <: IR_coeff}
    ir_is_valid_basic(ring, monoms, coeffs)
    for i in 1:length(monoms)
        !(length(monoms[i]) == length(coeffs[i])) && throw(DomainError("Invalid IR."))
        for j in 1:length(monoms[i])
            !(length(monoms[i][j]) == ring.nvars) && throw(DomainError("Invalid IR."))
            !(all(>=(0), monoms[i][j])) && throw(DomainError("Invalid IR."))
            iszero(coeffs[i][j]) && throw(DomainError("Invalid IR")) # can be relaxed
            if (ring.ground == :zp)
                !(0 < coeffs[i][j] < ring.characteristic) && throw(DomainError("Invalid IR."))
            end
        end
    end
    true
end

function ir_ensure_valid(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}}
) where {I <: Integer, C <: Coeff}
    ir_is_valid_basic(ring, monoms, coeffs)
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
            if ring.ground == :zp
                new_coeffs[i][j] = mod(new_coeffs[i][j], ring.characteristic)
            end
        end
    end
    # Sort terms if needed
    tdeg(e) = vcat(sum(e), e)
    for i in 1:length(new_monoms)
        if !issorted(new_monoms[i], lt=(a, b) -> monom_isless(tdeg(a), tdeg(b), ring.ord), rev=true)
            perm = collect(1:length(new_monoms[i]))
            sort!(
                perm,
                lt=(a, b) -> monom_isless(tdeg(new_monoms[i][a]), tdeg(new_monoms[i][b]), ring.ord),
                rev=true
            )
            new_monoms[i] = new_monoms[i][perm]
            new_coeffs[i] = new_coeffs[i][perm]
        end
    end
    # Merge terms
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
            if ring.ground == :zp && _new_coeffs[i][slow_idx] >= ring.characteristic
                _new_coeffs[i][slow_idx] -= ring.characteristic
                @invariant _new_coeffs[i][slow_idx] < ring.characteristic
            end
        end
    end
    new_monoms, new_coeffs = _new_monoms, _new_coeffs
    # Remove zero coefficients (zero polynomials stay)
    for i in 1:length(new_monoms)
        perm = collect(1:length(new_monoms[i]))
        filter!(j -> !iszero(new_coeffs[i][j]), perm)
        new_monoms[i] = new_monoms[i][perm]
        new_coeffs[i] = new_coeffs[i][perm]
    end
    ring, new_monoms, new_coeffs
end

# Process polynomials on the apply stage

function ir_extract_coeffs_raw!(trace, coeffs::Vector{Vector{C}}) where {C <: Coeff}
    # write new coefficients directly to trace buffer
    input_polys_perm = trace.input_permutation
    term_perms = trace.term_sorting_permutations
    homog_term_perms = trace.term_homogenizing_permutations
    CoeffsType = trace.representation.coefftype

    permute_input_terms = !isempty(term_perms)
    permute_homogenizing_terms = !isempty(homog_term_perms)

    @inbounds for i in 1:length(coeffs)
        basis_cfs = trace.buf_basis.coeffs[i]
        poly_index = input_polys_perm[i]
        poly = coeffs[poly_index]
        isempty(poly) && return false
        !(length(poly) == length(basis_cfs)) && return false
        for j in 1:length(poly)
            coeff_index = j
            if permute_input_terms
                coeff_index = term_perms[poly_index][coeff_index]
            end
            if permute_homogenizing_terms
                coeff_index = homog_term_perms[poly_index][coeff_index]
            end
            coeff = poly[coeff_index]
            basis_cfs[j] = convert(CoeffsType, coeff)
        end
    end

    if trace.params.homogenize
        @invariant trace.ring.ground == :zp
        trace.buf_basis.coeffs[length(coeffs) + 1][1] = one(CoeffsType)
        trace.buf_basis.coeffs[length(coeffs) + 1][2] =
            trace.ring.characteristic - one(typeof(trace.ring.characteristic))
    end

    true
end

# Packed coefficients utils

function ir_pack_coeffs(batch::NTuple{N, T}) where {N, T}
    ring = batch[1][1]
    ch = CompositeNumber(map(el -> el[1].characteristic, batch))
    new_ring = PolyRing(ring.nvars, ring.ord, ch, ring.ground)
    monoms = batch[1][2]
    coeffs = Vector{Vector{CompositeNumber{N, UInt64}}}(undef, length(monoms))
    @assert allequal(map(el -> el[2], batch))
    for i in 1:length(batch[1][2])
        coeffs[i] = Vector{CompositeNumber{N, UInt64}}(undef, length(batch[1][2][i]))
        for j in 1:length(batch[1][2][i])
            coeffs[i][j] = CompositeNumber(ntuple(k -> batch[k][3][i][j], N))
        end
    end
    true, new_ring, monoms, coeffs
end

function ir_unpack_composite_coefficients(
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

###
# Converting to internal representation

@timeit _TIMER function ir_convert_ir_to_internal(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    repr = params.representation
    monoms2 = Vector{Vector{repr.monomtype}}(undef, length(monoms))
    CT = repr.coefftype
    if !isconcretetype(CT)
        CT = C
    end

    coeffs2 = Vector{Vector{CT}}(undef, length(monoms))
    @inbounds for i in 1:length(monoms)
        monoms2[i] = Vector{repr.monomtype}(undef, length(monoms[i]))
        coeffs2[i] = Vector{CT}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            monoms2[i][j] = monom_construct_from_vector(repr.monomtype, monoms[i][j])
            coeffs2[i][j] = CT(coeffs[i][j])
        end
    end
    ring2, term_sorting_permutations = ir_set_monomial_ordering!(ring, monoms2, coeffs2, params)
    term_sorting_permutations, ring2, monoms2, coeffs2
end

@timeit _TIMER function ir_convert_internal_to_ir(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    monoms2 = Vector{Vector{Vector{IR_exponent}}}(undef, length(monoms))
    if eltype(eltype(coeffs)) <: AbstractFloat
        coeffs2 = map(cc -> map(c -> UInt(c), cc), coeffs)
    else
        coeffs2 = coeffs
    end
    @inbounds for i in 1:length(monoms)
        monoms2[i] = Vector{Vector{IR_exponent}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            monoms2[i][j] = Vector{IR_exponent}(undef, ring.nvars)
            monom_to_vector!(monoms2[i][j], monoms[i][j])
        end
    end
    monoms2, coeffs2
end

# Checks that the monomial ordering is consistent.
# Sorts the polynomials terms w.r.t. the target ordering.
function ir_set_monomial_ordering!(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params
) where {M <: Monom, C <: Coeff}
    ordering_check_consistency(ring.nvars, params.target_ord)
    if ring.ord == params.target_ord
        # No reordering of terms needed
        return ring, Vector{Vector{Int}}()
    end
    ring = PolyRing(ring.nvars, params.target_ord, ring.characteristic, ring.ground)
    permutations = sort_input_terms_to_change_ordering!(monoms, coeffs, params.target_ord)
    ring, permutations
end
