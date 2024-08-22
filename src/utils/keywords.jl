# This file is a part of Groebner.jl. License is GNU GPL v2.

# Handling keyword arguments for functions in the interface.

const INT_INF = typemax(Int)

#! format: off
# Syntax formatting is disabled for the next several lines.
#
# Maps a function name from the interface to a set of supported keyword
# arguments with their corresponding default values.
const _supported_kw_args = (
    groebner = (
        reduced      = true,
        ordering     = InputOrdering(),
        certify      = false,
        linalg       = :auto,
        monoms       = :auto,
        arithmetic   = :auto,
        seed         = 42,
        loglevel     = _loglevel_default,
        maxpairs     = INT_INF,
        selection    = :auto,
        modular      = :auto,
        threaded     = :auto,
        sweep        = false,
        homogenize   = :auto,
        statistics   = :no,
        batched      = true,
        changematrix = false
    ),
    normalform = (
        check       = false,
        ordering    = InputOrdering(),
        monoms      = :dense,
        loglevel    = _loglevel_default,
        statistics  = :no
    ),
    isgroebner = (
        ordering    = InputOrdering(),
        certify     = true,
        seed        = 42,
        monoms      = :dense,
        loglevel    = _loglevel_default,
        statistics  = :no
    ),
    groebner_learn = (
        seed        = 42,
        ordering    = InputOrdering(),
        monoms      = :auto,
        arithmetic  = :auto,
        loglevel    = _loglevel_default,
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
        loglevel    = _loglevel_default,
        sweep       = true,
        statistics  = :no,
        threaded    = :auto,
    ),
    groebner_with_change_matrix = (
        reduced      = true,
        ordering     = InputOrdering(),
        certify      = false,
        linalg       = :auto,
        monoms       = :auto,
        arithmetic   = :auto,
        seed         = 42,
        loglevel     = _loglevel_default,
        maxpairs     = INT_INF,
        selection    = :auto,
        modular      = :auto,
        threaded     = :auto,
        sweep        = false,
        homogenize   = :auto,
        statistics   = :no,
        batched      = true,
        changematrix = true
    ),
)
#! format: on

"""
    KeywordArguments

Stores keyword arguments passed to a function in the interface in Groebner.jl.
    
Can be manually created with 

```julia
kwargs = Groebner.KeywordArguments(:groebner, seed = 99, reduced = false)
```
"""
mutable struct KeywordArguments
    function_id::Symbol
    reduced::Bool
    ordering::Any
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
    changematrix::Bool

    KeywordArguments(function_id::Symbol; passthrough...) =
        KeywordArguments(function_id, passthrough)

    function KeywordArguments(function_id::Symbol, kws)
        @assert haskey(_supported_kw_args, function_id)
        default_kw_args = _supported_kw_args[function_id]
        for (key, _) in kws
            if !haskey(default_kw_args, key)
                error_msg = join(sort(map(string, collect(keys(default_kw_args)))), ", ")
                throw(
                    AssertionError("""
              Keyword \"$key\" is not supported by Groebner.$(function_id).
              Supported keyword arguments for Groebner.$(function_id) are:\n$error_msg""")
                )
            end
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
        @assert arithmetic in (:auto, :delayed, :signed, :basic, :floating) """
        Not recognized arithmetic: $arithmetic
        Possible choices for keyword "arithmetic" are:
        `:auto`, `:delayed`, `:signed`, `:basic`, `:floating`"""

        seed = get(kws, :seed, get(default_kw_args, :seed, 42))

        loglevel_int_or_sym =
            get(kws, :loglevel, get(default_kw_args, :loglevel, _loglevel_default))
        @assert loglevel_int_or_sym isa Integer ||
                loglevel_int_or_sym in _loglevels_spelled_out """
        Not recognized logging level: $loglevel_int_or_sym.
        Value passed to keyword "loglevel" must be either an integer or one of $_loglevels_spelled_out."""

        loglevel = if loglevel_int_or_sym in _loglevels_spelled_out
            _loglevel_spell_to_int[loglevel_int_or_sym]
        else
            loglevel_int_or_sym
        end

        maxpairs = get(kws, :maxpairs, get(default_kw_args, :maxpairs, INT_INF))
        @assert maxpairs > 0 "The limit on the number of critical pairs must be positive"

        modular = get(kws, :modular, get(default_kw_args, :modular, :auto))
        @assert modular in (:auto, :classic_modular, :learn_and_apply) """
        Not recognized modular strategy: $modular
        Possible choices for keyword "modular" are:
        `:auto`, `:classic_modular`, `:learn_and_apply`"""

        batched = get(kws, :batched, get(default_kw_args, :batched, true))

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

        changematrix = get(kws, :changematrix, get(default_kw_args, :changematrix, false))

        new(
            function_id,
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
            changematrix
        )
    end
end

function logging_setup(keywords::KeywordArguments)
    logger_update(loglevel=keywords.loglevel)
    nothing
end

function statistics_setup(keywords::KeywordArguments)
    log_simdinfo()
    nothing
end
