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

    # Use multi-threading.
    # This does nothing currently.
    threading::Bool

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
    #
    heuristic_check = true
    randomized_check = true
    certify_check = kwargs.certify
    # 
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
    #
    linalg = kwargs.linalg
    # Default linear algebra algorithm is randomized
    if linalg === :auto
        linalg = :randomized
    end
    if !iszero(ring.ch) && linalg === :randomized
        # Do not use randomized linear algebra if the field characteristic is
        # too small. 
        # TODO: In the future, it would be good to adapt randomized linear
        # algebra to this case by taking more random samples
        if ring.ch < 5000
            @log level = -1 """
            The field characteristic is too small.
            Switching from randomized linear algebra to a deterministic one."""
            linalg = :deterministic
        end
    end
    linalg_sparsity = :sparse
    linalg_algorithm = LinearAlgebra(linalg, linalg_sparsity)

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
    modular_strategy = kwargs.modular
    majority_threshold = 1
    #
    seed = kwargs.seed
    rng = _default_rng_type(seed)
    useed = UInt64(seed)
    #
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
    arithmetic = $arithmetic
    reduced = $reduced
    homogenize = $homogenize
    maxpairs = $maxpairs
    selection_strategy = $selection_strategy
    ground = $ground
    modular_strategy = $modular_strategy
    majority_threshold = $majority_threshold
    threading = $threading
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
        reduced,
        maxpairs,
        selection_strategy,
        ground,
        modular_strategy,
        majority_threshold,
        threading,
        useed,
        rng,
        sweep,
        statistics
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
        params.modular_strategy,
        params.majority_threshold,
        params.threading,
        params.seed,
        params.rng,
        params.sweep,
        params.statistics
    )
end
