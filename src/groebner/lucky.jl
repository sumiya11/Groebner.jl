# Lucky primes are prime numbers used in modular computation.
# The sequence of lucky primes starts with the FIRST_LUCKY_PRIME

# Good primes are prime numbers used for checking correctness
# (nice terminology)
# The sequence of good primes starts with the FIRST_GOOD_PRIME

# Groebner basis is computed modulo a lucky prime, and the correctness
# of the computed basis is checked modulo a good prime.

# There are 50697537 primes between FIRST_LUCKY_PRIME and FIRST_GOOD_PRIME.
# This limits us to about 50697537*log2(2^30)*log10(2) decimal digits of coefficient size.
#
# There are about twice as much prime numbers between FIRST_LUCKY_PRIME and 2^31-1.
# This limits us to about 2*50697537*log2(2^30)*log10(2) decimal digits of coefficient size.
#
# Hence, we are able to safely handle coefficients of up to 10^8 digits in size.
const FIRST_LUCKY_PRIME = 2^31 - 1
const FIRST_GOOD_PRIME  = 2^30 + 3

@noinline function _too_large_coefficient(modulo)
    throw(
        ErrorException(
            "Too large coefficient encountered in the basis (last prime is $modulo). This should not happen normally, please sumbit an issue."
        )
    )
end

# Keeps track of used prime numbers and helps selecting new ones
mutable struct LuckyPrimes
    # TODO in the future: once parametric coefficients are supported, revise this
    # integer coefficients of input polynomials
    coeffs::Vector{Vector{BigInt}}
    # buffer for BigInts
    buf::BigInt

    # current lucky prime
    luckyprime::UInt64
    # current good prime
    goodprime::UInt64

    # all used lucky primes
    primes::Vector{UInt64}
    # product of all used lucky primes
    modulo::BigInt

    function LuckyPrimes(coeffs::Vector{Vector{BigInt}})
        new(coeffs, BigInt(), FIRST_LUCKY_PRIME, FIRST_GOOD_PRIME, UInt64[], BigInt(1))
    end
end

updatemodulo!(tracker::LuckyPrimes) =
    Base.GMP.MPZ.mul_ui!(tracker.modulo, last(tracker.primes))

# Check if prime is lucky w.r.t. coefficients from tracker
# i.e., does not divide any of the leading coefficients
function isluckyprime(tracker::LuckyPrimes, prime::UInt64)
    buf = tracker.buf
    p   = BigInt(prime)
    for poly in tracker.coeffs
        for c in poly
            if Base.GMP.MPZ.cmp_ui(Base.GMP.MPZ.tdiv_r!(buf, c, p), 0) == 0
                return false
            end
        end
    end
    return true
end

# Returns the next lucky prime
function nextluckyprime!(tracker::LuckyPrimes)
    prime = tracker.luckyprime

    while !isluckyprime(tracker, prime)
        prime = nextprime(prime + 1)
    end

    if prime >= 2^32
        _too_large_coefficient(prime)
    end

    tracker.luckyprime = nextprime(prime + 1)
    push!(tracker.primes, prime)

    return prime
end

function nextgoodprime!(tracker::LuckyPrimes)
    prime = tracker.goodprime

    while !isluckyprime(tracker, prime)
        prime = nextprime(prime + 1)
    end

    if prime >= FIRST_LUCKY_PRIME
        _too_large_coefficient(prime)
    end

    tracker.goodprime = nextprime(prime + 1)

    return prime
end
