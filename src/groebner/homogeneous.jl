
function max_total_degree(ring::PolyRing, monoms::Vector{Vector{T}}) where {T}
    D = zero(T)
    @inbounds for i in 1:length(monoms)
        d = sum(monoms[i])
        D = max(d, D)
    end
    D
end

# function are_generators_homogeneous(
#     ring::PolyRing,
#     monoms::Vector{Vector{Vector{T}}},
#     coeffs::Vector{Vector{C}}
# ) where {T, C <: Coeff}
#     @assert length(monoms) == length(coeffs)
#     D = max_total_degree(ring, monoms)
#     @inbounds for i in 1:length(monoms)
#         for j in 1:length(monoms[i])
#             d = sum(monoms[i][j])
#             if d != D
#                 return false
#             end
#         end
#     end
#     true
# end

function homogenize_generators(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}},
    params
) where {T, C <: Coeff}
    @log level = -2 "Homogenizing generators"
    @assert length(monoms) == length(coeffs)
    nvars = ring.nvars
    new_nvars = nvars + 1
    new_monoms = Vector{Vector{Vector{T}}}(undef, length(monoms))
    @inbounds for i in 1:length(monoms)
        D = max_total_degree(ring, monoms[i])
        new_monoms[i] = Vector{Vector{T}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            new_monoms[i][j] = Vector{T}(undef, new_nvars)
            d = sum(monoms[i][j])
            @invariant d <= D
            new_monoms[i][j][end] = D - d
            for k in 1:nvars
                new_monoms[i][j][k] = monoms[i][j]
            end
        end
    end
    new_ord = _ProductOrdering(ring.ord, _Lex(Int[new_nvars]))
    new_ring = PolyRing(new_nvars, new_ord, ring.ch)
    @log level = -2 "Homogenized polynomial ring" new_ring
    new_ring, new_monoms, coeffs
end

function dehomogenize_generators(
    ring,
    monoms::Vector{Vector{Vector{T}}},
    coeffs::Vector{Vector{C}},
    params
) where {T, C <: Coeff}
    @assert length(monoms) == length(coeffs)
    @log level = -2 "De-homogenizing generators.."
    nvars = ring.nvars
    @assert nvars > 1
    new_nvars = nvars - 1
    new_monoms = Vector{Vector{Vector{T}}}(undef, length(monoms))
    @inbounds for i in 1:length(monoms)
        new_monoms[i] = Vector{Vector{T}}(undef, length(monoms[i]))
        for j in 1:length(monoms[i])
            new_monoms[i][j] = Vector{T}(undef, new_nvars)
            for k in 1:new_nvars
                new_monoms[i][j][k] = monoms[i][j]
            end
        end
    end
    new_ord = ring.ord.ord1
    new_ring = PolyRing(new_nvars, new_ord, ring.ch)
    @log level = -2 "De-homogenized polynomial ring" new_ring
    new_ring, new_monoms, coeffs
end

function saturate_generators(ring, monoms, coeffs, sat_monoms, sat_coeffs, params)
    1
end
