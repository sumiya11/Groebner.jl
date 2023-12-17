# Select parameters for Groebner basis computation

# It seems there is no Xoshiro rng in Julia v < 1.8.
# Use Random.Xoshiro, if available, as it is a bit faster.
const _default_rng_type = @static if VERSION >= v"1.8.0"
    Random.Xoshiro
else
    Random.MersenneTwister
end

# Specifies linear backend algorithm
struct LinearAlgebra
    # One of :deterministic, :randomized, :experimental_1, :experimental_2,
    # :experimental_3,
    algorithm::Symbol
    # One of :dense, :sparse
    sparsity::Symbol

    function LinearAlgebra(algorithm, sparsity)
        new(algorithm, sparsity)
    end
end

# Stores parameters for a single GB computation.
# NOTE: in principle, MonomOrd1, ..., MonomOrd3 can be subtypes of any type
# besides the usual Groebner.AbstractInternalOrdering
mutable struct AlgorithmParameters{
    MonomOrd1,
    MonomOrd2,
    MonomOrd3,
    Arithmetic <: AbstractArithmetic
}
    # Desired monomial ordering of output polynomials
    target_ord::MonomOrd1
    # Monomial ordering for the actual computation
    computation_ord::MonomOrd2
    # Original monomial ordering of input polynomials
    original_ord::MonomOrd3

    # Specifies correctness checks levels
    heuristic_check::Bool
    randomized_check::Bool
    certify_check::Bool

    # If do homogenize input generators
    homogenize::Bool

    # This option only makes sense for functions `normalform` and `kbase`. It
    # specifies if the program should check if the input is indeed a Groebner
    # basis.
    check::Bool

    # Linear algebra backend to be used
    linalg::LinearAlgebra

    # This can hold buffers or precomputed multiplicative inverses to speed up
    # the arithmetic in the ground field
    arithmetic::Arithmetic
    using_smallest_type_for_coeffs::Bool

    # If reduced Groebner basis is needed
    reduced::Bool

    # Limit the number of critical pairs in the F4 matrix by this number
    maxpairs::Int

    # Selection strategy. One of the following:
    # - :normal
    # well, it is tricky to implement sugar selection with F4..
    selection_strategy::Symbol

    # Ground field of computation. This can be one of the following:
    # - :qq for the rationals
    # - :zp for integers modulo a prime
    ground::Symbol

    # Strategy for modular computation in groebner. This can be one of the
    # following:
    # - :classic_modular
    # - :learn_and_apply
    modular_strategy::Symbol

    # In modular computation of the basis, compute (at least!) this many bases
    # modulo different primes until a consensus in majority vote is reached
    majority_threshold::Int

    # :simultaneous or :incremental
    crt_algorithm::Symbol

    # Use multi-threading.
    # This does nothing currently.
    threaded::Symbol

    # Random number generator
    seed::UInt64
    rng::_default_rng_type

    # Internal option for `groebner`.
    # At the end of F4, polynomials are interreduced. 
    # We can mark and sweep polynomials that are redundant prior to
    # interreduction to speed things up a bit. This option specifies if such
    # sweep should be done.
    sweep::Bool

    statistics::Symbol
end

function AlgorithmParameters(
    ring,
    representation,
    kwargs::KeywordsHandler;
    orderings=nothing
)
    # TODO: we should probably document this better
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

    # The levels of correctness checks. By default, we always check correctness
    # modulo a "random" prime
    heuristic_check = true
    randomized_check = true
    certify_check = kwargs.certify

    homogenize = if kwargs.homogenize === :yes
        true
    else
        if kwargs.homogenize === :auto
            if ring.nvars <= 1
                false
            elseif computation_ord isa Lex || computation_ord isa ProductOrdering
                true
            else
                false
            end
        else
            false
        end
    end

    linalg = kwargs.linalg
    if !iszero(ring.ch) && linalg === :randomized
        # Do not use randomized linear algebra if the field characteristic is
        # too small. 
        # TODO: In the future, it would be good to adapt randomized linear
        # algebra to this case by taking more random samples
        if ring.ch < 500
            @log level = 1_000 """
            The field characteristic is too small ($(ring.ch)).
            Switching from randomized linear algebra to a deterministic one."""
            linalg = :deterministic
        end
    end
    if linalg === :auto
        # Default linear algebra algorithm is randomized
        linalg = :randomized
    end
    # Default linear algebra algorithm is sparse
    linalg_sparsity = :sparse
    linalg_algorithm = LinearAlgebra(linalg, linalg_sparsity)

    arithmetic = select_arithmetic(
        ring.ch,
        representation.coefftype,
        kwargs.arithmetic,
        representation.using_smallest_type_for_coeffs
    )

    ground = :zp
    if iszero(ring.ch)
        ground = :qq
    end

    reduced = kwargs.reduced
    maxpairs = kwargs.maxpairs

    selection_strategy = kwargs.selection
    if selection_strategy === :auto
        if target_ord isa Union{Lex, ProductOrdering}
            selection_strategy = :normal # TODO :sugar
        else
            selection_strategy = :normal
        end
    end

    threaded = kwargs.threaded
    if !(_threaded[])
        if threaded === :yes
            @log level = 1_000 """
            You have explicitly provided the keyword argument `threaded = :yes`,
            however, the multi-threading is disabled globally in Groebner.jl due
            to the set environment variable GROEBNER_NO_THREADED=0

            Consider enabling threading by setting GROEBNER_NO_THREADED to 1"""
        end
        threaded = :no
    end

    # By default, modular computation uses learn & apply
    modular_strategy = kwargs.modular
    if modular_strategy === :auto
        modular_strategy = :learn_and_apply
    end
    if !reduced
        @log level = -1 """
        The option reduced=$reduced was passed in the input, 
        falling back to classic multi-modular algorithm."""
        modular_strategy = :classic_modular
    end

    majority_threshold = 1
    crt_algorithm = :simultaneous

    seed = kwargs.seed
    rng = _default_rng_type(seed)
    useed = UInt64(seed)

    sweep = kwargs.sweep

    statistics = kwargs.statistics

    @log level = -1 """
    Selected parameters:
    target_ord = $target_ord
    computation_ord = $computation_ord
    original_ord = $original_ord
    heuristic_check = $heuristic_check
    randomized_check = $randomized_check
    certify_check = $certify_check
    check = $(kwargs.check)
    linalg = $linalg_algorithm
    threaded = $threaded
    arithmetic = $arithmetic
    using_smallest_type_for_coeffs = $(representation.using_smallest_type_for_coeffs)
    reduced = $reduced
    homogenize = $homogenize
    maxpairs = $maxpairs
    selection_strategy = $selection_strategy
    ground = $ground
    modular_strategy = $modular_strategy
    majority_threshold = $majority_threshold
    crt_algorithm = $crt_algorithm
    seed = $seed
    rng = $rng
    sweep = $sweep
    statistics = $statistics"""

    AlgorithmParameters(
        target_ord,
        computation_ord,
        original_ord,
        heuristic_check,
        randomized_check,
        certify_check,
        homogenize,
        kwargs.check,
        linalg_algorithm,
        arithmetic,
        representation.using_smallest_type_for_coeffs,
        reduced,
        maxpairs,
        selection_strategy,
        ground,
        modular_strategy,
        majority_threshold,
        crt_algorithm,
        threaded,
        useed,
        rng,
        sweep,
        statistics
    )
end

function params_mod_p(
    params::AlgorithmParameters,
    prime::C;
    using_smallest_type_for_coeffs=nothing
) where {C <: Coeff}
    smallest_type_for_coeffs = if !isnothing(using_smallest_type_for_coeffs)
        using_smallest_type_for_coeffs
    else
        params.using_smallest_type_for_coeffs
    end
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
        select_arithmetic(prime, C, :auto, smallest_type_for_coeffs),
        smallest_type_for_coeffs,
        params.reduced,
        params.maxpairs,
        params.selection_strategy,
        params.ground,
        params.modular_strategy,
        params.majority_threshold,
        params.crt_algorithm,
        params.threaded,
        params.seed,
        params.rng,
        params.sweep,
        params.statistics
    )
end
