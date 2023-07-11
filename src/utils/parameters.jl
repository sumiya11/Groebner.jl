# Select parameters for Groebner basis computation

const _default_rng_type = @static if VERSION >= v"1.8.0"
    Random.Xoshiro
else
    Random.MersenneTwister
end

struct AlgorithmParameters{Ord1, Ord2}
    # Output polynomials monomial ordering
    target_ord::Ord1
    # Monomial ordering for computation
    computation_ord::Ord2

    # Basis correctness checks levels
    heuristic_check::Bool
    randomized_check::Bool
    certify_check::Bool

    check::Bool

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
    majority_threshold::Int
    emit_computation_graph::Bool

    threading::Bool

    # Random number generator seed
    seed::UInt64
    rng::_default_rng_type
end

function AlgorithmParameters(ring, kwargs::KeywordsHandler)
    ordering = kwargs.ordering
    if kwargs.ordering === InputOrdering() || kwargs.ordering === nothing
        ordering = ring.ord
    else
        ordering = kwargs.ordering
    end
    target_ord = ordering
    computation_ord = ordering

    heuristic_check = true
    randomized_check = true
    certify_check = kwargs.certify

    linalg = kwargs.linalg

    ground = :zp
    if iszero(ring.ch)
        ground = :qq
    end

    reduced = kwargs.reduced
    maxpairs = kwargs.maxpairs

    threading = false

    strategy = kwargs.strategy
    majority_threshold = 1

    seed = kwargs.seed

    emit_computation_graph = false

    rng = _default_rng_type(seed)

    useed = UInt64(seed)

    @log level = -1 """
    Selected parameters:
    target_ord = $target_ord
    computation_ord = $computation_ord
    heuristic_check = $heuristic_check
    randomized_check = $randomized_check
    certify_check = $certify_check
    check = $kwargs.check
    linalg = $linalg
    reduced = $reduced
    maxpairs = $maxpairs
    ground = $ground
    strategy = $strategy
    majority_threshold = $majority_threshold
    emit_computation_graph = $emit_computation_graph
    threading = $threading
    seed = $seed
    rng = $rng"""

    AlgorithmParameters(
        target_ord,
        computation_ord,
        heuristic_check,
        randomized_check,
        certify_check,
        kwargs.check,
        linalg,
        reduced,
        maxpairs,
        ground,
        strategy,
        majority_threshold,
        emit_computation_graph,
        threading,
        useed,
        rng
    )
end
