# Additional validation; the upstream SHA routine only retains its last basis polynomial.
import Nemo
include("../layman_parser.jl")

function full_basis_signature(text; normalize_monic=false)
    ring, basis = parse_system_naive(text)
    field = string(Nemo.characteristic(Nemo.base_ring(ring)))
    variables = string.(Nemo.gens(ring))
    # The parser constructs canonical Nemo polynomials in the same explicit
    # degrevlex order. Compare every coefficient and monomial, allowing only a
    # permutation of the output basis, not omitted elements. Some reference
    # backends emit primitive integer polynomials over QQ; explicitly opt in
    # to monic normalization when comparing those outputs.
    if normalize_monic
        basis = [iszero(f) ? f : Nemo.divexact(f, Nemo.leading_coefficient(f)) for f in basis]
    end
    return (field, variables, sort(string.(basis)))
end

function compare_full_bases(reference, candidate; normalize_monic=false)
    full_basis_signature(reference; normalize_monic) ==
        full_basis_signature(candidate; normalize_monic)
end

if abspath(PROGRAM_FILE) == @__FILE__
    normalize_monic = !isempty(ARGS) && first(ARGS) == "--normalize-monic"
    paths = normalize_monic ? ARGS[2:end] : ARGS
    length(paths) >= 2 || error("usage: validate.jl [--normalize-monic] REFERENCE_OUTPUT CANDIDATE_OUTPUT ...")
    reference = full_basis_signature(read(first(paths), String); normalize_monic)
    for candidate in paths[2:end]
        reference == full_basis_signature(read(candidate, String); normalize_monic) ||
            error("Full basis mismatch: $candidate")
        println("Full basis matches: $candidate")
    end
end
