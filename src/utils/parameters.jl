# Select parameters for groebner basis computation

struct AlgorithmParameters{Ord1, Ord2}
    # Output polynomials monomial ordering
    target_ord::Ord1
    # Monomial ordering for computation
    computation_ord::Ord2

    # Basis correctness checks levels
    heuristic_check::Bool
    randomized_check::Bool
    guaranteed_check::Bool

    # Linear algebra backend to be used. Currently available are
    # - :deterministic for exact deterministic algebra,
    # - :randomized for probabilistic linear algebra
    linalg::Symbol

    # Reduced Groebner basis is needed
    reduced::Bool

    maxpairs::Int

    # Ground field of computation. Currently options are
    # - :qq for rationals,
    # - :zp for integers modulo a prime.
    ground::Symbol

    # TODO
    strategy::Symbol

    threading::Bool

    # Random number generator seed
    seed::UInt64
end

function AlgorithmParameters(ring, kwargs::KeywordsHandler)
    ordering = kwargs.ordering
    if kwargs.ordering === InputOrdering() || kwargs.ordering === nothing
        ordering = ring.ord
    end
    target_ord = ordering
    computation_ord = ordering

    heuristic_check = true
    randomized_check = true
    guaranteed_check = kwargs.certify

    linalg = kwargs.linalg

    ground = :zp
    if iszero(ring.ch)
        ground = :qq
    end

    reduced = kwargs.reduced
    maxpairs = kwargs.maxpairs

    threading = false

    strategy = :classic_modular

    seed = kwargs.seed

    AlgorithmParameters(
        target_ord,
        computation_ord,
        heuristic_check,
        randomized_check,
        guaranteed_check,
        linalg,
        reduced,
        maxpairs,
        ground,
        strategy,
        threading,
        UInt64(seed)
    )
end
