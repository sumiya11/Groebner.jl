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
        linalg      = :randomized,
        monoms      = :auto,
        seed        = 42,
        loglevel    = _default_loglevel,
        maxpairs    = typemax(Int),   # NOTE: maybe use Inf?
        selection   = :auto,
        modular     = :classic_modular,
        sweep       = false,
        homogenize  = :auto,
        statistics  = :no
    ),
    normalform = (
        check    = false,
        ordering = InputOrdering(),
        monoms   = :dense,
        loglevel = _default_loglevel,
        statistics  = :no
    ),
    isgroebner = (
        ordering = InputOrdering(),
        certify  = true,
        seed     = 42,
        monoms   = :dense,
        loglevel = _default_loglevel,
        statistics  = :no
    ),
    kbase = (
        check    = false,
        ordering = InputOrdering(),
        monoms   = :dense,
        loglevel = _default_loglevel,
        statistics  = :no
    ),
    groebner_learn = (
        seed        = 42,
        ordering    = InputOrdering(),
        monoms      = :auto,
        loglevel    = _default_loglevel,
        homogenize  = :auto,
        sweep       = true,
        statistics  = :no
    ),
    groebner_apply! = (
        seed     = 42,
        monoms   = :auto,
        loglevel = _default_loglevel,
        sweep    = true,
        statistics  = :no
    )
)
#! format: on

"""
    KeywordsHandler

Stores keyword arguments passed to one of the functions in the interface. On
creation, checks that the arguments are correct.

Sets the global logger for Groebner.jl.
"""
struct KeywordsHandler{Ord}
    reduced::Bool
    ordering::Ord
    certify::Bool
    linalg::Symbol
    monoms::Symbol
    seed::Int
    loglevel::Int
    maxpairs::Int
    selection::Symbol
    modular::Symbol
    check::Bool
    sweep::Bool
    homogenize::Symbol
    statistics::Symbol

    function KeywordsHandler(function_key, kws)
        @assert haskey(_supported_kw_args, function_key)
        default_kw_args = _supported_kw_args[function_key]
        for (key, _) in kws
            @assert haskey(default_kw_args, key) "Not recognized keyword: $key"
        end
        reduced = get(kws, :reduced, get(default_kw_args, :reduced, true))
        ordering = get(kws, :ordering, get(default_kw_args, :ordering, InputOrdering()))
        certify = get(kws, :certify, get(default_kw_args, :certify, false))
        linalg = get(kws, :linalg, get(default_kw_args, :linalg, :randomized))
        @assert linalg in (:randomized, :deterministic) """
        Not recognized linear algebra option: $linalg
        Possible choices for keyword "linalg" are:
        `:randomized`, `:deterministic`"""
        monoms = get(kws, :monoms, get(default_kw_args, :monoms, :dense))
        @assert monoms in (:auto, :dense, :packed, :sparse) """
        Not recognized monomial representation: $monoms
        Possible choices for keyword "monoms" are:
        `:auto`, `:dense`, `:packed`, `:sparse`"""
        seed = get(kws, :seed, get(default_kw_args, :seed, 42))
        loglevel = get(kws, :loglevel, get(default_kw_args, :loglevel, 0))
        @assert loglevel isa Integer """
        Not recognized logging level: $loglevel::$(typeof(loglevel))
        Values passed to keyword "loglevel" must be integers."""
        update_logger(loglevel=loglevel)
        if loglevel <= 0
            refresh_performance_counters!()
        end
        maxpairs = get(kws, :maxpairs, get(default_kw_args, :maxpairs, typemax(Int)))
        @assert maxpairs > 0 "The limit on the number of critical pairs must be positive"
        modular = get(kws, :modular, get(default_kw_args, :modular, :classic_modular))
        @assert modular in (:classic_modular, :learn_and_apply) """
        Not recognized modular strategy: $modular
        Possible choices for keyword "modular" are:
        `:classic_modular`, `:learn_and_apply`"""
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
        @assert statistics in (:no, :timings, :all) """
        Not recognized option for collecting statistics: $statistics
        Possible choices for keyword "statistics" are:
        `:no`, `:timings`, `:all`"""
        @log level = -1 """
          Using keywords: 
          reduced    = $reduced, 
          ordering   = $ordering, 
          certify    = $certify, 
          linalg     = $linalg, 
          monoms     = $monoms, 
          seed       = $seed, 
          loglevel   = $loglevel, 
          maxpairs   = $maxpairs,
          selection  = $selection,
          modular   = $modular,
          check      = $check,
          sweep      = $sweep
          homogenize = $homogenize
          statistics = $statistics"""
        new{typeof(ordering)}(
            reduced,
            ordering,
            certify,
            linalg,
            monoms,
            seed,
            loglevel,
            maxpairs,
            selection,
            modular,
            check,
            sweep,
            homogenize,
            statistics
        )
    end
end
