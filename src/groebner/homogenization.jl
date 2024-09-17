# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Homogenization & saturation

function maximum_totaldeg(ring::PolyRing, monoms::Vector{Vector{T}}) where {T}
    D = zero(T)
    @inbounds for i in 1:length(monoms)
        d = monom_totaldeg(monoms[i])
        D = max(d, D)
    end
    D
end

function extend_ordering_in_homogenization(
    nvars::Int,
    ord::Ord,
    homogenizing_vairable=:last
) where {Ord <: Union{DegRevLex, DegLex, Lex, ProductOrdering}}
    @assert homogenizing_vairable === :last
    lex_part = Lex(nvars + 1)
    new_ord = ProductOrdering(ordering_make_not_simple(ord, nvars), lex_part)
    new_ord
end

function restrict_ordering_in_dehomogenization(ord, homogenizing_vairable=:last)
    ord.ord1
end

function extend_ordering_in_saturation(
    nvars::Int,
    ord::Ord,
    saturating_vairable=:last
) where {Ord <: Union{DegRevLex, DegLex, Lex, ProductOrdering}}
    @assert saturating_vairable === :last
    lex_part = Lex(nvars + 1)
    new_ord = ProductOrdering(lex_part, ordering_make_not_simple(ord, nvars))
    new_ord
end

function restrict_ordering_in_desaturation(ord, saturating_vairable=:last)
    ord.ord2
end

function homogenize_generators!(
    ring::PolyRing{Ord},
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}},
    params
) where {Ord, T, C <: Coeff}
    @assert length(monoms) == length(coeffs)
    nvars = ring.nvars
    new_nvars = nvars + 1
    new_monoms = Vector{Vector{Vector{T}}}(undef, length(monoms))
    @inbounds for i in 1:length(monoms)
        D = maximum_totaldeg(ring, monoms[i])
        new_monoms[i] = Vector{Vector{T}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            # `+ 1` since exponent vectors also store the total degree
            new_monoms[i][j] = Vector{T}(undef, new_nvars + 1)
            d = monom_totaldeg(monoms[i][j])
            @invariant d <= D
            new_monoms[i][j][1] = D
            new_monoms[i][j][end] = D - d
            for k in 2:(nvars + 1)
                new_monoms[i][j][k] = monoms[i][j][k]
            end
        end
    end
    # TODO: clarify the order of variables
    new_ord = extend_ordering_in_homogenization(ring.nvars, ring.ord)
    new_ring = PolyRing(new_nvars, new_ord, ring.ch)
    term_permutation = sort_input_terms_to_change_ordering!(new_monoms, coeffs, new_ord)
    sat_var_index = new_nvars
    new_ring_sat, new_monoms, coeffs =
        saturate_generators_by_variable!(new_ring, new_monoms, coeffs, params, sat_var_index)
    params.target_ord = new_ring_sat.ord
    term_permutation, new_ring_sat, new_monoms, coeffs
end

function dehomogenize_generators!(
    ring,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}},
    params
) where {T, C <: Coeff}
    ring_desat, monoms, coeffs = desaturate_generators!(ring, monoms, coeffs, params)
    @assert length(monoms) == length(coeffs)
    nvars = ring_desat.nvars
    @assert nvars > 1
    new_nvars = nvars - 1
    new_monoms = Vector{Vector{Vector{T}}}(undef, length(monoms))
    reduced_to_zero = Vector{Int}()
    @inbounds for i in 1:length(monoms)
        new_monoms[i] = Vector{Vector{T}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            new_monoms[i][j] = Vector{T}(undef, new_nvars + 1)
            for k in 2:(new_nvars + 1)
                new_monoms[i][j][k] = monoms[i][j][k]
            end
            new_monoms[i][j][1] = sum(new_monoms[i][j][2:end])
        end
        # if all(monom -> iszero(monom_totaldeg(monom)), new_monoms[i])
        #     push!(reduced_to_zero, i)
        # end
    end
    deleteat!(new_monoms, reduced_to_zero)
    deleteat!(coeffs, reduced_to_zero)
    new_ord = restrict_ordering_in_dehomogenization(ring_desat.ord)
    new_ring = PolyRing(new_nvars, new_ord, ring_desat.ch)
    sort_input_terms_to_change_ordering!(new_monoms, coeffs, new_ord)
    params.target_ord = new_ring.ord
    new_monoms, coeffs = autoreduce2(new_ring, new_monoms, coeffs, params)
    new_ring, new_monoms, coeffs
end

function desaturate_generators!(
    ring,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}},
    params
) where {T, C <: Coeff}
    @assert length(monoms) == length(coeffs)
    nvars = ring.nvars
    @assert nvars > 1
    new_nvars = nvars - 1
    new_monoms = Vector{Vector{Vector{T}}}(undef, length(monoms))
    new_sparse_row_coeffs = similar(coeffs)
    new_size = 0
    # remove the saturating variable
    @inbounds for i in 1:length(monoms)
        to_skip = false
        for j in 1:length(monoms[i])
            if monoms[i][j][end] > zero(T)
                to_skip = true
                break
            end
            to_skip && break
        end
        to_skip && break
        new_size += 1
        new_sparse_row_coeffs[new_size] = coeffs[i]
        new_monoms[new_size] = Vector{Vector{T}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            new_monoms[i][j] = Vector{T}(undef, new_nvars + 1)
            for k in 2:(new_nvars + 1)
                new_monoms[i][j][k] = monoms[i][j][k]
            end
            new_monoms[i][j][1] = sum(new_monoms[i][j][2:end])
        end
    end
    @assert new_size > 0
    resize!(new_sparse_row_coeffs, new_size)
    resize!(new_monoms, new_size)
    new_ord = restrict_ordering_in_desaturation(ring.ord)
    new_ring = PolyRing(new_nvars, new_ord, ring.ch)
    params.target_ord = new_ring.ord
    new_ring, new_monoms, new_sparse_row_coeffs
end

function saturate_generators_by_variable!(
    ring,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}},
    params,
    sat_var_index
) where {T, C <: Coeff}
    @assert length(monoms) == length(coeffs)
    nvars = ring.nvars
    new_nvars = nvars + 1
    new_monoms = Vector{Vector{Vector{T}}}(undef, length(monoms))
    @inbounds for i in 1:length(monoms)
        new_monoms[i] = Vector{Vector{T}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            # `+ 1` since exponent vectors also store the total degree
            new_monoms[i][j] = Vector{T}(undef, new_nvars + 1)
            for k in 1:(nvars + 1)
                new_monoms[i][j][k] = monoms[i][j][k]
            end
            new_monoms[i][j][end] = zero(T)
        end
    end
    # Construct `xi*t - 1`, where `i` is the `sat_var_index`
    new_poly_monoms = Vector{Vector{T}}(undef, 2)
    const_monom = zeros(T, new_nvars + 1)
    lead_monom = zeros(T, new_nvars + 1)
    lead_monom[sat_var_index + 1] = one(T)
    lead_monom[end] = one(T)
    lead_monom[1] = one(T) + one(T)
    new_poly_monoms[1] = lead_monom
    new_poly_monoms[2] = const_monom
    new_poly_coeffs = Vector{C}(undef, 2)
    new_poly_coeffs[1] = one(C)
    # NOTE: minus one in the current ground field
    new_poly_coeffs[2] = iszero(ring.ch) ? -one(C) : (ring.ch - one(ring.ch))
    push!(new_monoms, new_poly_monoms)
    push!(coeffs, new_poly_coeffs)
    new_ord = extend_ordering_in_saturation(ring.nvars, ring.ord)
    new_ring = PolyRing(new_nvars, new_ord, ring.ch)
    params.target_ord = new_ring.ord
    new_ring, new_monoms, coeffs
end
