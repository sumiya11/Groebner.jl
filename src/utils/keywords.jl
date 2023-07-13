# Handling keyword arguments in the interface.

const _default_loglevel = 0

#! format: off
# Syntax formatting is disabled for the next several lines.
#
# Maps a function name from the interface to a set of supported keyword
# arguments with their corresponding default values.
const _supported_kw_args = (
    groebner = (
        reduced  = true, 
        ordering = InputOrdering(),
        certify  = false,
        linalg   = :randomized,
        monoms   = :packed,
        seed     = 42,
        loglevel = _default_loglevel,
        maxpairs = typemax(Int),   # TODO,
        strategy = :classic_modular,
        sweep    = false,
    ),
    normalform = (
        check    = true,
        ordering = InputOrdering(),
        monoms   = :default,
        loglevel = _default_loglevel
    ),
    isgroebner = (
        ordering = InputOrdering(),
        certify  = true,
        seed     = 42,
        monoms   = :default,
        loglevel = _default_loglevel
    ),
    kbase = (
        check    = true,
        ordering = InputOrdering(),
        monoms   = :default,
        loglevel = _default_loglevel
    ),
    groebner_learn = (
        seed     = 42,
        monoms   = :packed,
        loglevel = _default_loglevel,
        sweep    = false,
    ),
    groebner_apply! = (
        seed     = 42,
        monoms   = :packed,
        loglevel = _default_loglevel,
        sweep    = false,
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
    strategy::Symbol
    check::Bool
    sweep::Bool

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
        @assert linalg in (:randomized, :deterministic) "Not recognized linear algebra option: $linalg"
        monoms = get(kws, :monoms, get(default_kw_args, :monoms, :packed))
        @assert monoms in (:default, :packed, :sparse) "Not recognized monomial representation: $monoms"
        seed = get(kws, :seed, get(default_kw_args, :seed, 42))
        loglevel = get(kws, :loglevel, get(default_kw_args, :loglevel, 0))
        update_logger(loglevel=loglevel)
        maxpairs = get(kws, :maxpairs, get(default_kw_args, :maxpairs, typemax(Int)))
        @assert maxpairs > 0 "The limit on the number of critical pairs must be positive"
        strategy = get(kws, :strategy, get(default_kw_args, :strategy, :classic_modular))
        @assert strategy in (:classic_modular, :learn_and_apply) "Not recognized strategy: $strategy"
        check = get(kws, :check, get(default_kw_args, :check, true))
        sweep = get(kws, :sweep, get(default_kw_args, :sweep, false))
        @log level = -1 """
          Using keywords: 
          reduced   = $reduced, 
          ordering  = $ordering, 
          certify   = $certify, 
          linalg    = $linalg, 
          monoms    = $monoms, 
          seed      = $seed, 
          loglevel  = $loglevel, 
          maxpairs  = $maxpairs,
          strategy  = $strategy,
          check     = $check,
          sweep     = $sweep"""
        new{typeof(ordering)}(
            reduced,
            ordering,
            certify,
            linalg,
            monoms,
            seed,
            loglevel,
            maxpairs,
            strategy,
            check,
            sweep
        )
    end
end
