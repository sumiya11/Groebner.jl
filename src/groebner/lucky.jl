# Lucky primes are prime numbers used in modular computation. 
#
# The sequence of lucky primes is increasing and deterministic, and starts with
# the FIRST_LUCKY_PRIME. Groebner basis is computed modulo a lucky prime, and
# the correctness of the computed basis is checked modulo some another prime.

# There are currently 50,697,537 primes between FIRST_LUCKY_PRIME and
# FIRST_CHECK_PRIME. This limits us to about 
#   50,697,537*log2(2^30)*log10(2) ~ 4 * 10^8 
# decimal digits of coefficient size. Hence, we are able to more or
# less safely handle coefficients of up to 10^8 digits.
const FIRST_LUCKY_PRIME = 2^31 - 1
const FIRST_CHECK_PRIME = 2^30 + 3

@noinline function __too_large_coefficient_error(modulo)
    throw(
        ErrorException(
            """
            Too large coefficient encountered in the basis (last used prime was $modulo). 
            This should not happen normally, please consider sumbitting an issue."""
        )
    )
end

# Keeps track of used prime numbers and helps selecting new ones
mutable struct LuckyPrimes
    # Integer coefficients of input polynomials
    coeffs::Vector{Vector{BigInt}}
    # buffer for BigInts
    buf::BigInt
    # current lucky prime
    luckyprime::UInt64
    # current check prime
    checkprime::UInt64
    # all used lucky primes
    primes::Vector{UInt64}
    # product of all used lucky primes
    modulo::BigInt

    function LuckyPrimes(coeffs::Vector{Vector{BigInt}})
        new(coeffs, BigInt(), FIRST_LUCKY_PRIME, FIRST_CHECK_PRIME, UInt64[], BigInt(1))
    end
end

# Check if the prime is lucky w.r.t. input basis coefficients -- does not divide
# any of the leading coefficients
function isluckyprime(lucky::LuckyPrimes, prime::UInt64)
    @log level = -3 "Checking if $prime is lucky.."
    buf = lucky.buf
    p   = BigInt(prime)
    for poly in lucky.coeffs
        for c in poly
            if Base.GMP.MPZ.cmp_ui(Base.GMP.MPZ.tdiv_r!(buf, c, p), 0) == 0
                return false
            end
        end
    end
    true
end

# Returns the next lucky prime
function next_lucky_prime!(lucky::LuckyPrimes)
    prime = lucky.luckyprime
    while !isluckyprime(lucky, prime)
        prime = Primes.nextprime(prime + 1)
        prime >= 2^32 && __too_large_coefficient_error(prime)
    end
    lucky.luckyprime = Primes.nextprime(prime + 1)
    @invariant !(prime in lucky.primes)
    push!(lucky.primes, prime)
    prime
end

function next_check_prime!(lucky::LuckyPrimes)
    prime = lucky.checkprime
    while !isluckyprime(lucky, prime)
        prime = nextprime(prime + 1)
        prime >= FIRST_LUCKY_PRIME && __too_large_coefficient_error(prime)
    end
    lucky.checkprime = nextprime(prime + 1)
    @invariant !(prime in lucky.primes)
    prime
end
