# Select parameters for Groebner basis computation

# It seems there is no Xoshiro rng in Julia v < 1.8.
# Use Random.Xoshiro, if available, as it is a bit faster.
const _default_rng_type = @static if VERSION >= v"1.8.0"
    Random.Xoshiro
else
    Random.MersenneTwister
end

# Stores parameters for a single GB computation.
mutable struct AlgorithmParameters{Ord1, Ord2, Ord3, Arithm}
    # Output polynomials monomial ordering
    target_ord::Ord1
    # Monomial ordering for computation
    computation_ord::Ord2
    # Original ordering
    original_ord::Ord3
    # Basis correctness checks levels
    heuristic_check::Bool
    randomized_check::Bool
    certify_check::Bool

    homogenize::Bool

    check::Bool

    # Linear algebra backend to be used. Currently available are
    # - :deterministic for exact deterministic algebra,
    # - :randomized for probabilistic linear algebra
    linalg::Symbol

    arithmetic::Arithm

    # Reduced Groebner basis is needed
    reduced::Bool

    maxpairs::Int

    selection_strategy::Symbol

    # Ground field of computation. Currently options are
    # - :qq for rationals,
    # - :zp for integers modulo a prime.
    ground::Symbol

    # TODO: introduce two strategies: :classic_modular and :learn_and_apply
    strategy::Symbol
    majority_threshold::Int

    threading::Bool

    # Random number generator seed
    seed::UInt64
    rng::_default_rng_type

    sweep::Bool
end

function AlgorithmParameters(
    ring,
    representation,
    kwargs::KeywordsHandler;
    orderings=nothing
)
    #
    if orderings !== nothing
        target_ord = orderings[2]
        computation_ord = orderings[2]
        original_ord = orderings[1]
    else
        if kwargs.ordering === InputOrdering() || kwargs.ordering === nothing
            ordering = ring.ord
        else
            ordering = kwargs.ordering
        end
        target_ord = ordering
        computation_ord = ordering
        original_ord = ring.ord
    end
    #
    heuristic_check = true
    randomized_check = true
    certify_check = kwargs.certify
    # 
    homogenize = if kwargs.homogenize === :yes
        true
    else
        if kwargs.homogenize === :auto
            computation_ord isa Lex || computation_ord isa ProductOrdering
        else
            false
        end
    end
    #
    linalg = kwargs.linalg
    arithmetic = select_arithmetic(ring.ch, representation.coefftype)
    ground = :zp
    if iszero(ring.ch)
        ground = :qq
    end
    #
    reduced = kwargs.reduced
    maxpairs = kwargs.maxpairs
    #
    selection_strategy = kwargs.selection
    if selection_strategy === :auto
        if target_ord isa Union{Lex, ProductOrdering}
            selection_strategy = :normal # :sugar
        else
            selection_strategy = :normal
        end
    end
    #
    threading = false
    #
    strategy = kwargs.strategy
    majority_threshold = 1
    #
    seed = kwargs.seed
    rng = _default_rng_type(seed)
    useed = UInt64(seed)
    #
    sweep = kwargs.sweep

    @log level = -1 """
    Selected parameters:
    target_ord = $target_ord
    computation_ord = $computation_ord
    original_ord = $original_ord
    heuristic_check = $heuristic_check
    randomized_check = $randomized_check
    certify_check = $certify_check
    check = $(kwargs.check)
    linalg = $linalg
    arithmetic = $arithmetic
    reduced = $reduced
    homogenize = $homogenize
    maxpairs = $maxpairs
    selection_strategy = $selection_strategy
    ground = $ground
    strategy = $strategy
    majority_threshold = $majority_threshold
    threading = $threading
    seed = $seed
    rng = $rng
    sweep = $sweep"""

    AlgorithmParameters(
        target_ord,
        computation_ord,
        original_ord,
        heuristic_check,
        randomized_check,
        certify_check,
        homogenize,
        kwargs.check,
        linalg,
        arithmetic,
        reduced,
        maxpairs,
        selection_strategy,
        ground,
        strategy,
        majority_threshold,
        threading,
        useed,
        rng,
        sweep
    )
end

function params_mod_p(params::AlgorithmParameters, prime::Integer)
    AlgorithmParameters(
        params.target_ord,
        params.computation_ord,
        params.original_ord,
        params.heuristic_check,
        params.randomized_check,
        params.certify_check,
        params.homogenize,
        params.check,
        params.linalg,
        select_arithmetic(prime, CoeffModular),
        params.reduced,
        params.maxpairs,
        params.selection_strategy,
        params.ground,
        params.strategy,
        params.majority_threshold,
        params.threading,
        params.seed,
        params.rng,
        params.sweep
    )
end
