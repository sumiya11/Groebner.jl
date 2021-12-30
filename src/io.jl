

struct PolyRing{Tv}
    #= Ring information =#
    # number of variables
    nvars::Int
    # raw length of exponent vector
    explen::Int
    # ring monomial ordering,
    # possible are :lex and :degrevlex
    ord::Symbol
end


# converts MPoly representation to internal polynomial representation used by algorithm
# by extracting base ring, exponents, and coefficients
#
function convert_to_internal(orig_polys::Vector{MPoly{Tv}}) where {Tv}
    # orig_polys is not empty here
    npolys = length(orig_polys)
    exps   = Vector{Vector{Vector{UInt16}}}(undef, npolys)
    coeffs = Vector{Vector{Tv}}(undef, npolys)

    for i in 1:npolys
        poly = orig_polys[i]
        exps[i] = Vector{Vector{UInt16}}(undef, length(poly))
        for j in 1:length(poly)
            # TODO: ask
            exps[i][j] = poly.exps[:, j]
        end
        coeffs[i] = copy(poly.coeffs)
    end

    R = parent(first(orig_polys))
    explen = R.N
    nvars  = R.num_vars
    ord    = R.ord

    ring = PolyRing{Tv}(nvars, explen, ord)

    println(exps, coeffs)
    return ring, exps, coeffs
end
