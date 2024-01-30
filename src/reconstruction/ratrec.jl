# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Single element rational reconstruction

"""
    ratrec_reconstruction_bound

Returns the bound for rational reconstruction based on the bitsize of the
modulo. As soon as the numerator in rational reconstruction exceeds this bound,
the gcd iteration stops.
"""
function ratrec_reconstruction_bound(modulo::BigInt)
    setprecision(2 * Base.GMP.MPZ.sizeinbase(modulo, 2)) do
        ceil(BigInt, sqrt(BigFloat(modulo) / 2))
    end
end

"""
    ratrec!

Computes the rational reconstruction of `a` mod `m`. Namely, a pair of numbers
`num`, `den`, such that 

    num//den ≡ a (mod m)

Writes the answer to `num` and `den` inplace. Returns `true` if the
reconstrction was successful and `false` otherwise. 

## Additional parameters:

- `bnd`: stores the stopping criterion threshold (see
    `ratrec_reconstruction_bound`) 
- `buf`, `buf1`, `buf2`, `buf3`, `u1`, `u2`, `u3`,  `v1`, `v2`, `v3`: buffers
"""
function ratrec!(
    num::BigInt,
    den::BigInt,
    bnd::BigInt,
    buf::BigInt,
    buf1::BigInt,
    buf2::BigInt,
    buf3::BigInt,
    u1::BigInt,
    u2::BigInt,
    u3::BigInt,
    v1::BigInt,
    v2::BigInt,
    v3::BigInt,
    a::BigInt,
    m::BigInt
)
    # Assumes the input is nonnegative!
    @invariant Base.GMP.MPZ.cmp_ui(a, 0) >= 0

    # fast path for numbers smaller than O(sqrt(modulo))
    if Base.GMP.MPZ.cmp(a, bnd) < 0
        Base.GMP.MPZ.set!(num, a)
        Base.GMP.MPZ.set_ui!(den, 1)
        return true
    end

    Base.GMP.MPZ.set_ui!(u1, 1)
    Base.GMP.MPZ.set_ui!(u2, 0)
    Base.GMP.MPZ.set!(u3, m)
    Base.GMP.MPZ.set_ui!(v1, 0)
    Base.GMP.MPZ.set_ui!(v2, 1)
    Base.GMP.MPZ.set!(v3, a)

    while true
        if Base.GMP.MPZ.cmp(v2, bnd) > 0
            return false
        end

        Base.GMP.MPZ.set!(buf, v3)
        if Base.GMP.MPZ.cmp_ui(buf, 0) < 0
            Base.GMP.MPZ.neg!(buf)
        end

        if Base.GMP.MPZ.cmp(buf, bnd) < 0
            break
        end

        Base.GMP.MPZ.tdiv_q!(buf, u3, v3)

        Base.GMP.MPZ.mul!(buf1, buf, v1)
        Base.GMP.MPZ.mul!(buf2, buf, v2)
        Base.GMP.MPZ.mul!(buf3, buf, v3)

        Base.GMP.MPZ.sub!(buf1, u1, buf1)
        Base.GMP.MPZ.sub!(buf2, u2, buf2)
        Base.GMP.MPZ.sub!(buf3, u3, buf3)

        Base.GMP.MPZ.set!(u1, v1)
        Base.GMP.MPZ.set!(u2, v2)
        Base.GMP.MPZ.set!(u3, v3)

        Base.GMP.MPZ.set!(v1, buf1)
        Base.GMP.MPZ.set!(v2, buf2)
        Base.GMP.MPZ.set!(v3, buf3)
    end

    Base.GMP.MPZ.set!(den, v2)
    Base.GMP.MPZ.set!(num, v3)

    #=
    Base.GMP.MPZ.gcd!(buf, den, num)
    Base.GMP.MPZ.tdiv_q!(den, buf)
    Base.GMP.MPZ.tdiv_q!(num, buf)
    =#

    if Base.GMP.MPZ.cmp_ui(den, 0) < 0
        Base.GMP.MPZ.neg!(den)
        Base.GMP.MPZ.neg!(num)
    end

    true
end

function ratrec_nemo(a::Nemo.ZZRingElem, m::Nemo.ZZRingElem)
    success, pq = Nemo.unsafe_reconstruct(a, m)
    success, (BigInt(numerator(pq)), BigInt(denominator(pq)))
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

    @inbounds for k in 1:length(indices)
        i, j = indices[k]
        rem_nemo = Nemo.ZZRingElem(table_zz[i][j])

        success, (num, den) = ratrec_nemo(rem_nemo, modulo_nemo)
        table_qq[i][j] = Base.unsafe_rational(num, den)

        !success && return false
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
    @assert length(table_qq) == length(table_zz)
    @assert modulo > 1

    modulo_nemo = Nemo.ZZRingElem(modulo)
    @inbounds for i in 1:length(table_zz)
        @assert length(table_zz[i]) == length(table_qq[i])
        for j in 1:length(table_zz[i])
            rem_nemo = Nemo.ZZRingElem(table_zz[i][j])
            @assert rem_nemo >= 0

            success, (num, den) = ratrec_nemo(rem_nemo, modulo_nemo)
            table_qq[i][j] = Base.unsafe_rational(num, den)

            !success && return false
        end
    end

    true
end
