


# given sparse vector of coefficients wrt monomial basis monom_basis
# constructs and returns corresponding polynomial (almost) inplace
function unsafe_sparsevector_to_poly(
                vector::SparseVectorAA{Tv, Ti},
                monom_basis) where {Tv, Ti}

    par  = parent(first(monom_basis))
    len  = nnz(vector)
    raw_explen = size(first(monom_basis).exps, 1)

    cfs  = Vector{Tv}(undef, len)
    exps = Matrix{UInt64}(undef, raw_explen, len)

    for (i, midx, val) in zip(1:len, nonzeroinds(vector), nonzeros(vector))
        cfs[i] = val
        for eidx in 1:raw_explen
            exps[eidx, i] = monom_basis[midx].exps[eidx, 1]
        end
    end

    MPoly{Tv}(par, cfs, exps)
end
