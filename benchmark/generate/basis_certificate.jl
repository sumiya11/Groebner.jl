import Nemo
using SHA

include((@__DIR__) * "/layman_parser.jl")

function hash_polynomial(poly, h)
    poly_str = repr(poly)
    sha3_256(poly_str)
end

function get_certificate(ring, polys)
    filter!(!iszero, polys)
    @assert !isempty(polys)
    sort!(polys, by=Nemo.leading_monomial)
    polys = map(poly -> Nemo.divexact(poly, Nemo.leading_coefficient(poly)), polys)
    h = map(i -> UInt8(0), 1:32)
    for (i, poly) in enumerate(polys)
        i_8 = i % UInt8
        h = i_8 .⊻ hash_polynomial(poly, h)
    end
    certificate = bytes2hex(h)
    certificate
end

function is_certificate_standardized(certificate)
    certificate == standardize_certificate(certificate)
end

function standardize_certificate(certificate)
    certificate = strip(string(certificate), ['\n', ' ', '\r'])
    certificate
end

function compute_basis_validation_certificate(result)
    ring, polys = parse_system_naive(result)
    if isempty(polys)
        return false, ""
    end
    certificate = get_certificate(ring, polys)
    true, standardize_certificate(certificate)
end
