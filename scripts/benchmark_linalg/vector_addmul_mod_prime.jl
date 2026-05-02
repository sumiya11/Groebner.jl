using AbstractAlgebra
using BitIntegers
using Nemo
using Primes
using Printf
using Random
using Statistics

const ROW_LENGTH = 50_000
const SUPPORT_LENGTH = 5_000
const PRIME64 = UInt64(18446744073709551557)
const PRIME256_BIG = Primes.prevprime(big(2)^255)
const PRIME256 = UInt256(PRIME256_BIG)
const PRIME512 = UInt512(PRIME256_BIG)
const PRIME25519_BIG = big(2)^255 - 19
const PRIME25519 = UInt256(PRIME25519_BIG)
const PRIME64_WIDE = UInt128(PRIME64)
const PRIME64_INV = Base.multiplicativeinverse(PRIME64_WIDE)
const PRIME256_INV = Base.multiplicativeinverse(PRIME512)
const PRIME25519_MASK = (one(UInt512) << 255) - one(UInt512)

function rand_big_lt(rng::AbstractRNG, p::BigInt)
    x = (BigInt(rand(rng, UInt128)) << 128) + BigInt(rand(rng, UInt128))
    mod(x, p)
end

function make_indices(rng::AbstractRNG, n::Int, k::Int)
    sort!(randperm(rng, n)[1:k])
end

@inline addmod_u64(a::UInt64, b::UInt64) = UInt64(rem(UInt128(a) + UInt128(b), PRIME64_INV))
@inline mulmod_u64(a::UInt64, b::UInt64) = UInt64(rem(widen(a) * widen(b), PRIME64_INV))

@inline addmulmod_u256(a::UInt256, mul::UInt256, coeff::UInt256) =
    UInt256(rem(UInt512(a) + UInt512(mul) * UInt512(coeff), PRIME256_INV))

@inline function reduce_curve25519(x::UInt512)
    y = (x & PRIME25519_MASK) + UInt512(19) * (x >> 255)
    z = (y & PRIME25519_MASK) + UInt512(19) * (y >> 255)
    r = UInt256(z)
    if r >= PRIME25519
        r -= PRIME25519
    end
    r
end

@inline addmulmod_curve25519(a::UInt256, mul::UInt256, coeff::UInt256) =
    reduce_curve25519(UInt512(a) + UInt512(mul) * UInt512(coeff))

function groebner_like_addmul_u64!(row::Vector{UInt64}, indices::Vector{Int}, coeffs::Vector{UInt64})
    mul = PRIME64 - row[indices[1]]
    @inbounds for j in eachindex(indices, coeffs)
        idx = indices[j]
        row[idx] = addmod_u64(row[idx], mulmod_u64(mul, coeffs[j]))
    end
    row
end

function groebner_like_addmul_u256!(row::Vector{UInt256}, indices::Vector{Int}, coeffs::Vector{UInt256})
    mul = PRIME256 - row[indices[1]]
    @inbounds for j in eachindex(indices, coeffs)
        idx = indices[j]
        row[idx] = addmulmod_u256(row[idx], mul, coeffs[j])
    end
    row
end

function groebner_like_addmul_curve25519!(
    row::Vector{UInt256},
    indices::Vector{Int},
    coeffs::Vector{UInt256}
)
    mul = PRIME25519 - row[indices[1]]
    @inbounds for j in eachindex(indices, coeffs)
        idx = indices[j]
        row[idx] = addmulmod_curve25519(row[idx], mul, coeffs[j])
    end
    row
end

function groebner_like_addmul_bigint!(row::Vector{BigInt}, indices::Vector{Int}, coeffs::Vector{BigInt}, p::BigInt)
    mul = p - row[indices[1]]
    @inbounds for j in eachindex(indices, coeffs)
        idx = indices[j]
        row[idx] = mod(row[idx] + mul * coeffs[j], p)
    end
    row
end

function groebner_like_addmul_generic!(row, indices::Vector{Int}, coeffs)
    mul = -row[indices[1]]
    @inbounds for j in eachindex(indices, coeffs)
        idx = indices[j]
        row[idx] = row[idx] + mul * coeffs[j]
    end
    row
end

function build_inputs(; row_length::Int=ROW_LENGTH, support_length::Int=SUPPORT_LENGTH, seed::Int=0)
    rng = MersenneTwister(seed)
    indices = make_indices(rng, row_length, support_length)

    row64 = [rand(rng, UInt64) % PRIME64 for _ in 1:row_length]
    coeffs64 = [rand(rng, UInt64) % PRIME64 for _ in 1:support_length]
    coeffs64[1] = one(UInt64)

    row_big = [rand_big_lt(rng, PRIME256_BIG) for _ in 1:row_length]
    coeffs_big = [rand_big_lt(rng, PRIME256_BIG) for _ in 1:support_length]
    coeffs_big[1] = BigInt(1)

    row_25519_big = [rand_big_lt(rng, PRIME25519_BIG) for _ in 1:row_length]
    coeffs_25519_big = [rand_big_lt(rng, PRIME25519_BIG) for _ in 1:support_length]
    coeffs_25519_big[1] = BigInt(1)

    row_u256 = UInt256.(row_big)
    coeffs_u256 = UInt256.(coeffs_big)
    row_25519_u256 = UInt256.(row_25519_big)
    coeffs_25519_u256 = UInt256.(coeffs_25519_big)

    F256 = Nemo.GF(PRIME256_BIG)
    row_nemo = F256.(row_big)
    coeffs_nemo = F256.(coeffs_big)
    coeffs_nemo[1] = one(F256)

    R256, _ = AbstractAlgebra.residue_ring(AbstractAlgebra.ZZ, PRIME256_BIG)
    row_aa = R256.(row_big)
    coeffs_aa = R256.(coeffs_big)
    coeffs_aa[1] = one(R256)

    F64 = Nemo.GF(BigInt(PRIME64))
    row_nemo64 = F64.(BigInt.(row64))
    coeffs_nemo64 = F64.(BigInt.(coeffs64))
    coeffs_nemo64[1] = one(F64)

    (
        ;
        indices,
        row64,
        coeffs64,
        row_nemo64,
        coeffs_nemo64,
        row_big,
        coeffs_big,
        row_u256,
        coeffs_u256,
        row_25519_big,
        coeffs_25519_big,
        row_25519_u256,
        coeffs_25519_u256,
        row_nemo,
        coeffs_nemo,
        row_aa,
        coeffs_aa,
    )
end

function verify_correctness(data)
    touched = data.indices

    row_big = copy(data.row_big)
    groebner_like_addmul_bigint!(row_big, data.indices, data.coeffs_big, PRIME256_BIG)
    @assert iszero(row_big[first(touched)])

    row_u256 = copy(data.row_u256)
    groebner_like_addmul_u256!(row_u256, data.indices, data.coeffs_u256)
    @assert iszero(row_u256[first(touched)])
    @assert all(BigInt(row_u256[idx]) == row_big[idx] for idx in touched)

    row_25519_big = copy(data.row_25519_big)
    groebner_like_addmul_bigint!(row_25519_big, data.indices, data.coeffs_25519_big, PRIME25519_BIG)
    @assert iszero(row_25519_big[first(touched)])

    row_25519_u256 = copy(data.row_25519_u256)
    groebner_like_addmul_curve25519!(row_25519_u256, data.indices, data.coeffs_25519_u256)
    @assert iszero(row_25519_u256[first(touched)])
    @assert all(BigInt(row_25519_u256[idx]) == row_25519_big[idx] for idx in touched)

    row_nemo = copy(data.row_nemo)
    groebner_like_addmul_generic!(row_nemo, data.indices, data.coeffs_nemo)
    @assert iszero(row_nemo[first(touched)])
    @assert all(BigInt(lift(Nemo.ZZ, row_nemo[idx])) == row_big[idx] for idx in touched)

    row_aa = copy(data.row_aa)
    groebner_like_addmul_generic!(row_aa, data.indices, data.coeffs_aa)
    @assert iszero(row_aa[first(touched)])

    row64 = copy(data.row64)
    groebner_like_addmul_u64!(row64, data.indices, data.coeffs64)
    @assert iszero(row64[first(touched)])

    row_nemo64 = copy(data.row_nemo64)
    groebner_like_addmul_generic!(row_nemo64, data.indices, data.coeffs_nemo64)
    @assert iszero(row_nemo64[first(touched)])
end

function run_case(
    name::AbstractString,
    kernel!,
    row_seed,
    updates_per_call;
    warmup::Int=5,
    samples::Int=30
)
    for _ in 1:warmup
        row = copy(row_seed)
        kernel!(row)
    end

    GC.gc()
    timings_ns = Vector{Float64}(undef, samples)
    for i in 1:samples
        row = copy(row_seed)
        t0 = time_ns()
        kernel!(row)
        timings_ns[i] = time_ns() - t0
    end

    median_ns = median(timings_ns)
    min_ns = minimum(timings_ns)
    updates_per_second = 1.0e9 * updates_per_call / median_ns
    (; name, median_ns, min_ns, updates_per_second)
end

function main()
    data = build_inputs()
    verify_correctness(data)

    results = [
        run_case(
            "UInt64 + MultInv",
            row -> groebner_like_addmul_u64!(row, data.indices, data.coeffs64),
            data.row64,
            length(data.indices)
        ),
        run_case(
            "Nemo GF(64-bit)",
            row -> groebner_like_addmul_generic!(row, data.indices, data.coeffs_nemo64),
            data.row_nemo64,
            length(data.indices)
        ),
        run_case(
            "BigInt mod p",
            row -> groebner_like_addmul_bigint!(row, data.indices, data.coeffs_big, PRIME256_BIG),
            data.row_big,
            length(data.indices)
        ),
        run_case(
            "UInt256 + MultInv",
            row -> groebner_like_addmul_u256!(row, data.indices, data.coeffs_u256),
            data.row_u256,
            length(data.indices)
        ),
        run_case(
            "UInt256 + 2^255-19",
            row -> groebner_like_addmul_curve25519!(row, data.indices, data.coeffs_25519_u256),
            data.row_25519_u256,
            length(data.indices)
        ),
        run_case(
            "Nemo GF(256-bit)",
            row -> groebner_like_addmul_generic!(row, data.indices, data.coeffs_nemo),
            data.row_nemo,
            length(data.indices)
        ),
        run_case(
            "AA residue(256-bit)",
            row -> groebner_like_addmul_generic!(row, data.indices, data.coeffs_aa),
            data.row_aa,
            length(data.indices)
        ),
    ]

    println("Groebner-like sparse/dense addmul benchmark")
    println("row_length=$(ROW_LENGTH), support_length=$(SUPPORT_LENGTH)")
    println(
        "prime64_bits=$(ndigits(BigInt(PRIME64), base=2)), generic_u256_prime_bits=$(ndigits(PRIME256_BIG, base=2)), special_prime=2^255-19"
    )
    println()
    @printf("%-22s %14s %14s %16s\n", "backend", "median ns", "best ns", "updates/sec")
    for result in results
        @printf(
            "%-22s %14.0f %14.0f %16.2f\n",
            result.name,
            result.median_ns,
            result.min_ns,
            result.updates_per_second
        )
    end
end

main()
