# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Intermediate representation of polynomials (ir)

# Polynomials in the intermediate representation are represented by two arrays:
# an array of exponent vectors and an array of coefficients. Zero polynomial is
# represented with two empty arrays.

# Intermediate representation (ir) enforces a number assumtpions:
# - Intermediate representation owns the memory.
# - Monomials in intermediate representation are sorted and unique.
# - Coefficients in intermediate representation are nonzero and normalized.

# Our conventions on some edge cases:
# - Trying to compute a Groebner basis of an empty set is an error.
# - The Groebner basis of [0] is [0].
# - The Groebner basis of [f1,...,fn, 0] is the Groebner basis of [f1...fn]

"""
    PolyRing

A polynomial ring.

## Example

4 variables, integers modulo 2^30 + 3, degrevlex order.

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

Base.:(==)(r1::PolyRing, r2::PolyRing) = r1.nvars == r2.nvars && r1.ord == r2.ord && r1.ch == r2.ch

ir_basic_is_valid(batch) = throw(DomainError("Invalid IR, unknown types."))
ir_basic_is_valid(ring, monoms, coeffs) = throw(DomainError("Invalid IR, unknown types."))

function ir_basic_is_valid(batch::NTuple{N, T}) where {N, T}
    for el in batch
        ir_basic_is_valid(el...)
    end
end

function ir_basic_is_valid(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}}
) where {T <: Integer, C <: Number}
    !(length(monoms) == length(coeffs)) && throw(DomainError("Invalid IR."))
    isempty(monoms) && throw(DomainError("Invalid IR."))
    !(ring.nvars >= 0) && throw(DomainError("The number of variables must be non-negative."))
    !(ring.ch >= 0) && throw(DomainError("Field characteristic must be nonnegative."))
    if ring.ch > 0
        !(C <: Integer) && throw(DomainError("Coefficients must be integers."))
        (C <: BigInt) && throw(DomainError("Coefficients must fit in a machine register."))
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
    # Sort terms if needed
    tdeg(e) = vcat(sum(e), e)
    for i in 1:length(new_monoms)
        if !issorted(
            new_monoms[i],
            lt=(a, b) -> monom_isless(tdeg(a), tdeg(b), ring.ord),
            rev=true
        )
            perm = collect(1:length(new_monoms[i]))
            sort!(
                perm,
                lt=(a, b) ->
                    monom_isless(tdeg(new_monoms[i][a]), tdeg(new_monoms[i][b]), ring.ord),
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
            if ring.ch > 0 && _new_coeffs[i][slow_idx] >= ring.ch
                _new_coeffs[i][slow_idx] -= ring.ch
                @invariant _new_coeffs[i][slow_idx] < ring.ch
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

# Packed coefficients utils

function ir_pack_coeffs(batch::NTuple{N, T}) where {N, T}
    ring = batch[1][1]
    ch = CompositeNumber(map(el -> el[1].ch, batch))
    new_ring = PolyRing(ring.nvars, ring.ord, ch)
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

function ir_convert_ir_to_internal(ring, monoms, coeffs, params, repr)
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
        ir_set_monomial_ordering!(ring, monoms2, coeffs2, params)
    term_sorting_permutations, ring2, monoms2, coeffs2
end

function ir_convert_internal_to_ir(ring, monoms, coeffs, params)
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

# Checks that the monomial ordering is consistent.
# Sorts the polynomials terms w.r.t. the target ordering.
function ir_set_monomial_ordering!(ring, monoms, coeffs, params)
    ordering_check_consistency(ring.nvars, params.target_ord)
    if ring.ord == params.target_ord
        # No reordering of terms needed
        return ring, Vector{Vector{Int}}()
    end
    ring = PolyRing(ring.nvars, params.target_ord, ring.ch)
    permutations = sort_input_terms_to_change_ordering!(monoms, coeffs, params.target_ord)
    ring, permutations
end
