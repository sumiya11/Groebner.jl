
#------------------------------------------------------------------------------

#=
    Lucky primes are prime numbers used for modular computation.
    The sequence of lucky primes starts with FIRST_LUCKY_PRIME

    Good primes are prime numbers used for checking correctness.
    The sequence of good primes starts with FIRST_GOOD_PRIME

    There are 50697537 primes between
        FIRST_LUCKY_PRIME and FIRST_GOOD_PRIME
=#

const FIRST_LUCKY_PRIME = 2^31 - 1
const FIRST_GOOD_PRIME  = 2^30 + 3

mutable struct PrimeTracker
    # integer generator coefficients from the input
    coeffs::Vector{Vector{CoeffZZ}}
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

    function PrimeTracker(coeffs_zz::Vector{Vector{CoeffZZ}})
        new(coeffs_zz, BigInt(), FIRST_LUCKY_PRIME, FIRST_GOOD_PRIME,
            UInt64[], BigInt(1))
    end
end

#------------------------------------------------------------------------------

# check if `prime` is lucky w.r.t. coefficients from tracker
function isluckyprime(tracker::PrimeTracker, prime::UInt64)
    buf = tracker.buf
    p   = BigInt(prime)
    for poly in tracker.coeffs
        for c in poly
            if Base.GMP.MPZ.tdiv_r!(buf, c, p) == 0
                 return false
            end
        end
    end
    return true
end

function nextluckyprime!(tracker::PrimeTracker)
    prime = tracker.luckyprime

    while !isluckyprime(tracker, prime)
        prime = nextprime(prime + 1)
    end

    tracker.luckyprime = nextprime(prime + 1)
    push!(tracker.primes, prime)

    return prime
end

#------------------------------------------------------------------------------

function updatemodulo!(tracker::PrimeTracker)
    Base.GMP.MPZ.mul_ui!(tracker.modulo, last(tracker.primes))
end

#------------------------------------------------------------------------------

function nextgoodprime!(tracker::PrimeTracker)
    prime = tracker.goodprime

    while !isluckyprime(tracker, prime)
        prime = nextprime(prime + 1)
    end

    if prime >= FIRST_LUCKY_PRIME
        error("Too large coefficient encountered in the basis. This should not happen normally")
    end

    tracker.goodprime = nextprime(prime + 1)

    return prime
end

#------------------------------------------------------------------------------
