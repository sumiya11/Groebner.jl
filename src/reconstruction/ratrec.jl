# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Single element rational reconstruction

# Returns the largest integer N (possibly off by 1) such that 2 N^2 < m.
function ratrec_reconstruction_bound(m::BigInt)
    isqrt((m >> UInt64(1)) - 1)
end

function ratrec_nemo(a::Nemo.ZZRingElem, m::Nemo.ZZRingElem)
    success, pq = Nemo.unsafe_reconstruct(a, m)
    success, Rational{BigInt}(pq)
end

function ratrec_nemo_2(
    a::Nemo.ZZRingElem,
    m::Nemo.ZZRingElem,
    N::Nemo.ZZRingElem,
    D::Nemo.ZZRingElem
)
    success, pq = Nemo.reconstruct(a, m, N, D)
    success, Rational{BigInt}(pq)
end

###
# Element-wise rational reconstruction

"""
    ratrec_vec_partial!(table_qq, table_zz, modulo, indices)

Given `table_zz` that stores remainders w.r.t. big `modulo`, reconstructs a
table of rationals using rational reconstruction.

Writes the resulting table to `table_qq`.
Returns `true` if all reconstructions succeeded and `false`, otherwise.

Elements of `table_zz` must be non-negative.
!!! Elements of `table_qq` must be initialized.

## Arguments

- `table_qq`: a table (vector of vectors) of big rationals. The result will be
  written here.
- `tables_zz`: a table of big integers, remainders w.r.t. `modulo`.
- `modulo`: a big integer.
- `indices`: a vector of tuples (int, int). Each tuple is an index to the table
  of remainders. If `(i, j)` is present in `indices`, then `table_qq[i][j]` will
  be reconstructed. Otherwise, the value of `table_qq[i][j]` is untouched.

## Example

Consider the following example

```julia
# Buffer
table_qq = [[BigInt(0) // BigInt(1)], [BigInt(0) // BigInt(1), BigInt(0) // BigInt(1)]]

# Residuals
table_zz = [[BigInt(58)], [BigInt(15), BigInt(73)]]
modulo = BigInt(77)

# Indices in the table to be reconstructed
indices = [(1, 1), (2, 1)]
success = Groebner.ratrec_vec_partial!(table_qq, table_zz, modulo, indices)
```

As a result, we obtain

```julia
@show table_qq success;
# table_qq = Vector{Rational{BigInt}}[[1//4], [-2//5, 0]]
# success = true
```
"""
function ratrec_vec_partial!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    indices::Vector{Tuple{Int, Int}}
)
    @invariant length(table_qq) == length(table_zz)
    modulo_nemo = Nemo.ZZRingElem(modulo)
    bnd = ratrec_reconstruction_bound(modulo)
    nemo_bnd = Nemo.ZZRingElem(bnd)

    @inbounds for k in 1:length(indices)
        i, j = indices[k]
        rem_nemo = Nemo.ZZRingElem(table_zz[i][j])

        success, pq = ratrec_nemo_2(rem_nemo, modulo_nemo, nemo_bnd, nemo_bnd)
        !success && return false

        table_qq[i][j] = pq
    end

    true
end

"""
    ratrec_vec_full!

Same as `ratrec_vec_partial!`, but reconstructs for all indices.
"""
function ratrec_vec_full!(
    table_qq::Vector{Vector{Rational{BigInt}}},
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt
)
    # indices = [(j, i) for j in 1:length(table_zz) for i in 1:length(table_zz[j])]
    # ratrec_vec_partial!(table_qq, table_zz, modulo, indices)
    @invariant length(table_qq) == length(table_zz)
    @invariant modulo > 1
    bnd = ratrec_reconstruction_bound(modulo)
    nemo_bnd = Nemo.ZZRingElem(bnd)

    modulo_nemo = Nemo.ZZRingElem(modulo)
    @inbounds for i in 1:length(table_zz)
        @invariant length(table_zz[i]) == length(table_qq[i])
        for j in 1:length(table_zz[i])
            rem_nemo = Nemo.ZZRingElem(table_zz[i][j])
            @invariant 0 <= rem_nemo < modulo

            success, pq = ratrec_nemo_2(rem_nemo, modulo_nemo, nemo_bnd, nemo_bnd)
            !success && return false

            table_qq[i][j] = pq
        end
    end

    true
end
