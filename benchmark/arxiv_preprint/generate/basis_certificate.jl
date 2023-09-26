include("myeval.jl")

import Nemo

function parse_system(result::String)
    lines = split(result, "\n")
    @assert length(lines) > 2
    vars_str = map(strip, split(lines[1], ","))
    char = parse(BigInt, strip(lines[2]))
    @assert char < typemax(UInt)
    polys_str = map(s -> strip(s, [' ', ',']), lines[3:end])
    polys_str = filter(!isempty, polys_str)
    base_field = iszero(char) ? Nemo.QQ : Nemo.GF(UInt(char))
    ring_nemo, vars_nemo = Nemo.PolynomialRing(base_field, vars_str, ordering=:degrevlex)
    polys_nemo = Vector{elem_type(ring_nemo)}()
    var_mapping = Dict{Symbol, elem_type(ring_nemo)}(
        Symbol(x) => gen(ring_nemo, i) for (i, x) in enumerate(vars_str)
    )
    for poly_str in polys_str
        poly_expr = Meta.parse(poly_str)
        poly_nemo = myeval(poly_expr, base_field, var_mapping)
        push!(polys_nemo, poly_nemo)
    end
    ring_nemo, polys_nemo
end

function get_certificate(ring, polys)
    length(polys)
end

function is_certificate_standardized(certificate)
    certificate == standardize_certificate(certificate)
end

function standardize_certificate(certificate)
    certificate = strip(string(certificate), ['\n', ' '])
    certificate
end

function compute_basis_validation_hash(result)
    ring, polys = parse_system(result)
    certificate = get_certificate(ring, polys)
    standardize_certificate(certificate)
end
