# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Lucky primes 

# The sequence of lucky prime candidates is decreasing and deterministic.
# Groebner basis is computed modulo a lucky prime, and the correctness of the
# computed basis is checked modulo some another prime(s).
#
# We use 31 bit primes for modular computation, since it allows the use of
# signed integers for implementing arithmetic in the ground field.

const FIRST_LUCKY_PRIME_CANDIDATE = 2^31 - 1
const FIRST_AUX_PRIME_CANDIDATE = 2^30 + 3

# Keeps track of used prime numbers and helps selecting new ones
mutable struct LuckyPrimes
    # Integer coefficients of input polynomials
    coeffs::Vector{Vector{BigInt}}
    buf::BigInt
    lucky_prime::UInt64
    aux_prime::UInt64
    used_primes::Vector{UInt64}
    # product of used lucky primes
    modulo::BigInt

    function LuckyPrimes(coeffs::Vector{Vector{BigInt}})
        new(
            coeffs,
            BigInt(),
            FIRST_LUCKY_PRIME_CANDIDATE,
            FIRST_AUX_PRIME_CANDIDATE,
            UInt64[],
            BigInt(1)
        )
    end
end

# Check if the prime is lucky, that is, if the prime does not divide any of the
# leading coefficients of the input.
function primes_is_lucky_prime(lucky::LuckyPrimes, prime::UInt64)
    buf = lucky.buf
    p   = BigInt(prime)
    @inbounds for coeffs in lucky.coeffs
        @invariant !isempty(coeffs)
        # TODO(Sasha): I am too greedy to check all coefficients.
        # for c in coeffs
        #     if Base.GMP.MPZ.cmp_ui(Base.GMP.MPZ.tdiv_r!(buf, c, p), 0) == 0
        #         return false
        #     end
        # end
        c = coeffs[1]
        if Base.GMP.MPZ.cmp_ui(Base.GMP.MPZ.tdiv_r!(buf, c, p), 0) == 0
            return false
        end
        c = coeffs[end]
        if Base.GMP.MPZ.cmp_ui(Base.GMP.MPZ.tdiv_r!(buf, c, p), 0) == 0
            return false
        end
    end
    true
end

function primes_next_lucky_prime!(lucky::LuckyPrimes)
    prime = lucky.lucky_prime
    while !primes_is_lucky_prime(lucky, prime)
        prime = Primes.prevprime(prime - 1)
    end
    lucky.lucky_prime = Primes.prevprime(prime - 1)
    @invariant !(prime in lucky.used_primes)
    push!(lucky.used_primes, prime)
    prime
end

function primes_next_aux_prime!(lucky::LuckyPrimes)
    prime = lucky.aux_prime
    while !primes_is_lucky_prime(lucky, prime)
        prime = nextprime(prime + 1)
        prime >= FIRST_LUCKY_PRIME_CANDIDATE && throw(
            DomainError(
                "Something went wrong in Groebner.jl. Please consider submitting a GitHub issue."
            )
        )
    end
    lucky.aux_prime = nextprime(prime + 1)
    @invariant !(prime in lucky.used_primes)
    prime
end
