# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Single element CRT

"""
    crt!

Implements Chinese Remainder lifting. Computes the unique `x` such that
        
    x ≡ a1 mod m1,
    x ≡ a2 mod m2.

Writes the answer to `buf` inplace.

## Additional parameters:

- `M`: must be equal to `m1 * m2`
- `buf`, `n1`, `n2`: additional buffers
- `c1`, `c2`: numbers, such that `c1 = m2 * invmod(m2, m1)`, and `c2 = m1 *
  invmod(m1, m2)`

Then, `x` is obtained as `x = c1 a1 + c2 a2 mod M`.
"""
function crt!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    a1::BigInt,
    c1::BigInt,
    a2::UInt,
    c2::BigInt
)
    Base.GMP.MPZ.mul_ui!(n1, c2, a2)
    Base.GMP.MPZ.mul!(n2, c1, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

"""
    crt_precompute!

Given the two moduli `m1` and `m2`, precomputes the multipliers `c1`, `c2` and
the modulo `M` for CRT.

`n1` and `n2` are additional buffers.
"""
function crt_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    c1::BigInt,
    m1::BigInt,
    c2::BigInt,
    m2::UInt
)
    Base.GMP.MPZ.mul_ui!(M, m1, m2)
    Base.GMP.MPZ.set_ui!(n2, m2)

    Base.GMP.MPZ.gcdext!(n1, c1, c2, n2, m1)
    Base.GMP.MPZ.mul_ui!(c1, m2)
    Base.GMP.MPZ.mul!(c2, m1)

    nothing
end

# Same as crt!, but for the case when `a2` is large
function crt!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    a1::BigInt,
    c1::BigInt,
    a2::BigInt,
    c2::BigInt
)
    Base.GMP.MPZ.mul!(n1, c2, a2)
    Base.GMP.MPZ.mul!(n2, c1, a1)

    Base.GMP.MPZ.add!(buf, n1, n2)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

# Same as crt_precompute!, but for the case when `m2` is large
function crt_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    c1::BigInt,
    m1::BigInt,
    c2::BigInt,
    m2::BigInt
)
    Base.GMP.MPZ.mul!(M, m1, m2)

    Base.GMP.MPZ.gcdext!(n1, c1, c2, m2, m1)
    Base.GMP.MPZ.mul!(c1, m2)
    Base.GMP.MPZ.mul!(c2, m1)

    nothing
end

"""
    crt!

Implements Chinese Remainder lifting with multiple moduli. Computes the unique
`x` such that for all `i` we have
        
    x ≡ ai mod mi.

Writes the answer to `buf` inplace.

## Additional parameters:

- `M`: must be equal to `∏ mi`
- `buf`, `n1`, `n2`: additional buffers
- `ci`: an array of numbers, with `ci[i] = πi * invmod(πi, mi)`, `πi = M / mi`

Then, `x` is obtained as `x = ∑ ci[i] ai[i] mod M`.
"""
function crt!(
    M::BigInt,
    buf::BigInt,
    n1::BigInt,
    n2::BigInt,
    ai::Vector{UInt},
    ci::Vector{BigInt}
)
    @invariant length(ai) == length(ci)

    Base.GMP.MPZ.set_ui!(n1, UInt(0))
    for i in 1:length(ai)
        Base.GMP.MPZ.mul_ui!(n2, ci[i], ai[i])
        Base.GMP.MPZ.add!(n1, n2)
    end

    Base.GMP.MPZ.set!(buf, n1)
    Base.GMP.MPZ.fdiv_r!(buf, M)

    nothing
end

"""
    crt_precompute!

Given the moduli, precomputes the multipliers `ci`, and the modulo `M` for CRT.

`n1` and `n2` are additional buffers.
"""
function crt_precompute!(
    M::BigInt,
    n1::BigInt,
    n2::BigInt,
    ci::Vector{BigInt},
    moduli::Vector{UInt}
)
    @invariant length(ci) == length(moduli)

    n3, n4 = BigInt(), BigInt()
    @inbounds Base.GMP.MPZ.set_ui!(M, moduli[1])
    @inbounds for i in 2:length(moduli)
        Base.GMP.MPZ.mul_ui!(M, moduli[i])
    end

    @inbounds for i in 1:length(moduli)
        Base.GMP.MPZ.set_ui!(n2, moduli[i])
        Base.GMP.MPZ.tdiv_q!(ci[i], M, n2)
        Base.GMP.MPZ.gcdext!(n2, n3, n4, ci[i], n2)
        Base.GMP.MPZ.mul!(ci[i], n3)
    end

    nothing
end

###
# Element-wise CRT

"""
    crt_vec_partial!(table_zz, modulo, tables_ff, moduli, indices)

Given `tables_ff` that stores remainders w.r.t. `moduli`, reconstructs a single
table of big integers using CRT.

Writes the resulting table to `table_zz` and the resulting big modulo to
`modulo`.

Each of `moduli` must fit in `UInt64`.
Elements of `tables_ff` must be non-negative.
!!! Elements of `table_zz` must be initialized and disjoint in memory (for
example, each element can be initilized by `BigInt(0)`).

## Arguments

- `table_zz`: a table (vector of vectors) of big integers. The result will be
  written here.
- `modulo`: a big integer. This function will modify `modulo` so that `modulo =
  prod(moduli)`.
- `tables_ff`: a vector of tables of machine integers. `tables_ff[i]` must be a
  table that stores remainders modulo `moduli[i]`.
- `moduli`: a vector of machine integers.
- `indices`: a vector of tuples (int, int). Each tuple is an index to the table
  of remainders. If `(i, j)` is present in `indices`, then `table_zz[i][j]` will
  be reconstructed. Otherwise, the value of `table_zz[i][j]` is untouched.

## Example

Consider the following example

```julia
# Buffers
table_zz = [[BigInt(0)], [BigInt(0), BigInt(0)]]
modulo = BigInt(0)

# Remainders
tables_ff = [
    [UInt64[2], UInt64[1, 3]],
    [UInt64[3], UInt64[4, 7]],
]
moduli = UInt64[7, 11]

# Indices in the table to be reconstructed
indices = [(1, 1), (2, 1)]

Groebner.crt_vec_partial!(table_zz, modulo, tables_ff, moduli, indices)
```

As a result, we obtain

```
@show table_zz modulo;
# table_zz = Vector{BigInt}[[58], [15, 0]]
# modulo = 77
```

Indeed, `77 = 7 * 11`, and 
    
    58 mod 7 == 2,
    58 mod 11 == 3.

Moreover, the element at index (2, 2) has not been reconstructed.
"""
function crt_vec_partial!(
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    tables_ff::Vector{Vector{Vector{T}}},
    moduli::Vector{T},
    indices::Vector{Tuple{Int, Int}}
) where {T <: Integer}
    @invariant isbitstype(T)
    n = length(moduli)
    @invariant n > 0

    # Base case
    if n == 1
        table_ff = tables_ff[1]
        Base.GMP.MPZ.set_ui!(modulo, UInt64(moduli[1]))
        @invariant length(table_zz) == length(table_ff)
        @inbounds for k in 1:length(indices)
            i, j = indices[k]
            rem_ij = UInt64(table_ff[i][j])
            Base.GMP.MPZ.set_ui!(table_zz[i][j], rem_ij)
        end
        return nothing
    end

    # n > 1
    # Precompute CRT multipliers
    buf = BigInt()
    n1, n2 = BigInt(), BigInt()
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    crt_precompute!(modulo, n1, n2, mults, map(UInt64, moduli))

    # Reconstruct an integer for each index (i, j) 
    rems = Vector{UInt64}(undef, n)
    @inbounds for k in 1:length(indices)
        i, j = indices[k]

        for t in 1:n
            rems[t] = UInt64(tables_ff[t][i][j])
        end

        crt!(modulo, buf, n1, n2, rems, mults)

        Base.GMP.MPZ.set!(table_zz[i][j], buf)
    end

    nothing
end

"""
    crt_vec_full!

Same as `crt_vec_partial!`, but reconstructs for all indices.
"""
function crt_vec_full!(
    table_zz::Vector{Vector{BigInt}},
    modulo::BigInt,
    tables_ff::Vector{Vector{Vector{T}}},
    moduli::Vector{T}
) where {T <: Integer}
    # indices = [(j, i) for j in 1:length(table_zz) for i in 1:length(table_zz[j])]
    # crt_vec_partial!(table_zz, modulo, tables_ff, moduli, indices)
    @assert isbitstype(T)
    @assert length(moduli) > 0

    n = length(moduli)
    # Base case
    if n == 1
        table_ff = tables_ff[1]
        Base.GMP.MPZ.set_ui!(modulo, UInt64(moduli[1]))
        @inbounds for i in 1:length(table_zz)
            @assert length(table_zz[i]) == length(table_ff[i])
            for j in 1:length(table_zz[i])
                rem_ij = UInt64(table_ff[i][j])
                @assert 0 <= rem_ij < moduli[1]
                table_zz[i][j] = rem_ij
            end
        end
        return nothing
    end

    # n > 1
    # Precompute CRT multipliers
    buf = BigInt()
    n1, n2 = BigInt(), BigInt()
    mults = Vector{BigInt}(undef, n)
    for i in 1:length(mults)
        mults[i] = BigInt(0)
    end
    crt_precompute!(modulo, n1, n2, mults, map(UInt64, moduli))

    rems = Vector{UInt64}(undef, n)
    @inbounds for i in 1:length(table_zz)
        for j in 1:length(table_zz[i])
            for t in 1:n
                @assert length(table_zz[i]) == length(tables_ff[t][i])
                @assert 0 <= tables_ff[t][i][j] < moduli[t]
                rems[t] = UInt64(tables_ff[t][i][j])
            end

            crt!(modulo, buf, n1, n2, rems, mults)

            Base.GMP.MPZ.set!(table_zz[i][j], buf)
        end
    end

    nothing
end
