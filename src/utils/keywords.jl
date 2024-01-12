# This file is a part of Groebner.jl. License is GNU GPL v2.

# Handling keyword arguments in the interface

const _default_loglevel = 0

#! format: off
# Syntax formatting is disabled for the next several lines.
#
# Maps a function name from the interface to a set of supported keyword
# arguments with their corresponding default values.
const _supported_kw_args = (
    groebner = (
        reduced     = true,
        ordering    = InputOrdering(),
        certify     = false,
        linalg      = :auto,
        monoms      = :auto,
        arithmetic  = :auto,
        seed        = 42,
        loglevel    = _default_loglevel,
        maxpairs    = typemax(Int),   # NOTE: maybe use Inf?
        selection   = :auto,
        modular     = :auto,
        threaded    = :auto,
        sweep       = false,
        homogenize  = :auto,
        statistics  = :no,
        batched     = true,
        use_flint   = true
    ),
    normalform = (
        check       = false,
        ordering    = InputOrdering(),
        monoms      = :dense,
        loglevel    = _default_loglevel,
        statistics  = :no
    ),
    isgroebner = (
        ordering    = InputOrdering(),
        certify     = true,
        seed        = 42,
        monoms      = :dense,
        loglevel    = _default_loglevel,
        statistics  = :no
    ),
    kbase = (
        check       = false,
        ordering    = InputOrdering(),
        monoms      = :dense,
        loglevel    = _default_loglevel,
        statistics  = :no
    ),
    fglm = (
        check       = false,
        monoms      = :dense,
        loglevel    = _default_loglevel,
    ),
    groebner_learn = (
        seed        = 42,
        ordering    = InputOrdering(),
        monoms      = :auto,
        arithmetic  = :auto,
        loglevel    = _default_loglevel,
        homogenize  = :auto,
        sweep       = true,
        statistics  = :no,
        threaded    = :auto,
    ),
    groebner_apply! = (
        seed        = 42,
        ordering    = InputOrdering(),
        monoms      = :auto,
        arithmetic  = :auto,
        loglevel    = _default_loglevel,
        sweep       = true,
        statistics  = :no,
        threaded    = :auto,
    )
)
#! format: on

"""
    KeywordsHandler

Stores keyword arguments passed to one of the functions in the interface."""
struct KeywordsHandler{Ord}
    reduced::Bool
    ordering::Ord
    certify::Bool
    linalg::Symbol
    threaded::Symbol
    monoms::Symbol
    arithmetic::Symbol
    seed::Int
    loglevel::Int
    maxpairs::Int
    selection::Symbol
    modular::Symbol
    batched::Bool
    check::Bool
    sweep::Bool
    homogenize::Symbol
    statistics::Symbol
    use_flint::Bool

    function KeywordsHandler(function_key, kws)
        @assert haskey(_supported_kw_args, function_key)
        default_kw_args = _supported_kw_args[function_key]
        for (key, _) in kws
            @assert haskey(default_kw_args, key) "Not recognized keyword: $key"
        end

        reduced = get(kws, :reduced, get(default_kw_args, :reduced, true))
        ordering = get(kws, :ordering, get(default_kw_args, :ordering, InputOrdering()))
        certify = get(kws, :certify, get(default_kw_args, :certify, false))

        linalg = get(kws, :linalg, get(default_kw_args, :linalg, :auto))
        @assert linalg in (
            :auto,
            :randomized,
            :deterministic,
            :experimental_1,
            :experimental_2,
            :experimental_3
        ) """
        Not recognized linear algebra option: $linalg
        Possible choices for keyword "linalg" are:
        `:auto`, `:randomized`, `:deterministic`"""

        threaded = get(kws, :threaded, get(default_kw_args, :threaded, :auto))
        @assert threaded in (:auto, :no, :yes, :force_yes) """
        Not recognized threading option: $threaded
        Possible choices for keyword "threaded" are:
        `:auto`, `:no`, `:yes`"""

        monoms = get(kws, :monoms, get(default_kw_args, :monoms, :dense))
        @assert monoms in (:auto, :dense, :packed, :sparse) """
        Not recognized monomial representation: $monoms
        Possible choices for keyword "monoms" are:
        `:auto`, `:dense`, `:packed`, `:sparse`"""

        arithmetic = get(kws, :arithmetic, get(default_kw_args, :arithmetic, :auto))
        @assert arithmetic in (:auto, :delayed, :signed, :basic) """
        Not recognized arithmetic: $arithmetic
        Possible choices for keyword "arithmetic" are:
        `:auto`, `:delayed`, `:signed`, `:basic`"""

        seed = get(kws, :seed, get(default_kw_args, :seed, 42))
        loglevel = get(kws, :loglevel, get(default_kw_args, :loglevel, 0))
        @assert loglevel isa Integer """
        Not recognized logging level: $loglevel::$(typeof(loglevel))
        Values passed to keyword "loglevel" must be integers."""

        maxpairs = get(kws, :maxpairs, get(default_kw_args, :maxpairs, typemax(Int)))
        @assert maxpairs > 0 "The limit on the number of critical pairs must be positive"

        modular = get(kws, :modular, get(default_kw_args, :modular, :auto))
        @assert modular in (:auto, :classic_modular, :learn_and_apply) """
        Not recognized modular strategy: $modular
        Possible choices for keyword "modular" are:
        `:auto`, `:classic_modular`, `:learn_and_apply`"""

        batched = get(kws, :batched, get(default_kw_args, :batched, true))
        use_flint = get(kws, :use_flint, get(default_kw_args, :use_flint, true))

        selection = get(kws, :selection, get(default_kw_args, :selection, :auto))
        @assert selection in (:auto, :normal, :sugar, :be_divided_and_perish)

        check = get(kws, :check, get(default_kw_args, :check, true))
        sweep = get(kws, :sweep, get(default_kw_args, :sweep, false))
        homogenize = get(kws, :homogenize, get(default_kw_args, :homogenize, :auto))
        @assert homogenize in (:auto, :no, :yes) """
        Not recognized homogenization strategy: $homogenize
        Possible choices for keyword "homogenize" are:
        `:auto`, `:no`, `:yes`"""

        statistics = get(kws, :statistics, get(default_kw_args, :statistics, :no))
        @assert statistics in (:no, :timings, :stats, :all) """
        Not recognized option for collecting statistics: $statistics
        Possible choices for keyword "statistics" are:
        `:no`, `:timings`, `:stats`, `:all`"""

        # @log level = -2 """
        #   In $function_key, using keywords: 
        #   reduced    = $reduced, 
        #   ordering   = $ordering, 
        #   certify    = $certify, 
        #   linalg     = $linalg, 
        #   monoms     = $monoms, 
        #   seed       = $seed, 
        #   loglevel   = $loglevel, 
        #   maxpairs   = $maxpairs,
        #   selection  = $selection,
        #   modular    = $modular,
        #   check      = $check,
        #   sweep      = $sweep
        #   homogenize = $homogenize
        #   statistics = $statistics"""

        new{typeof(ordering)}(
            reduced,
            ordering,
            certify,
            linalg,
            threaded,
            monoms,
            arithmetic,
            seed,
            loglevel,
            maxpairs,
            selection,
            modular,
            batched,
            check,
            sweep,
            homogenize,
            statistics,
            use_flint
        )
    end
end

function logging_setup(keywords::KeywordsHandler)
    logger_update(loglevel=keywords.loglevel)
    nothing
end

function statistics_setup(keywords::KeywordsHandler)
    if keywords.loglevel <= 0
        performance_counters_refresh()
        statistics_refresh()
    end
    nothing
end
