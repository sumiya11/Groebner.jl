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
        selection    = :auto,
        modular      = :auto,
        threaded     = :auto,
        homogenize   = :auto,
        changematrix = false,
        _composite   = 4,
        _generic     = false
    ),
    normalform = (
        check       = false,
        ordering    = InputOrdering(),
        monoms      = :dense,
    ),
    lead = (
        ordering    = InputOrdering(),
    ),
    isgroebner = (
        ordering    = InputOrdering(),
        certify     = true,
        seed        = 42,
        monoms      = :dense
    ),
    groebner_learn = (
        seed        = 42,
        ordering    = InputOrdering(),
        monoms      = :auto,
        arithmetic  = :auto,
        homogenize  = :auto,
        threaded    = :auto,
    ),
    groebner_apply! = (
        seed        = 42,
        ordering    = InputOrdering(),
        monoms      = :auto,
        arithmetic  = :auto,
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
        selection    = :auto,
        modular      = :auto,
        threaded     = :auto,
        homogenize   = :auto,
        changematrix = true,
        _composite   = 4,
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
    selection::Symbol
    modular::Symbol
    _composite::Int
    check::Bool
    homogenize::Symbol
    changematrix::Bool
    _generic::Bool

    function KeywordArguments(function_id::Symbol, kws)
        @assert haskey(_supported_kw_args, function_id)
        default_kw_args = _supported_kw_args[function_id]
        for (key, _) in kws
            if !haskey(default_kw_args, key)
                error_msg = join(sort(map(string, collect(keys(default_kw_args)))), ", ")
                throw(AssertionError("""
                Keyword \"$key\" is not supported by Groebner.$(function_id).
                Supported keyword arguments for Groebner.$(function_id) are:\n$error_msg"""))
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
        @assert monoms in (:auto, :dense, :packed) """
        Not recognized monomial representation: $monoms
        Possible choices for keyword "monoms" are:
        `:auto`, `:dense`, `:packed`"""

        arithmetic = get(kws, :arithmetic, get(default_kw_args, :arithmetic, :auto))
        @assert arithmetic in (:auto, :delayed, :signed, :basic) """
        Not recognized arithmetic: $arithmetic
        Possible choices for keyword "arithmetic" are:
        `:auto`, `:delayed`, `:signed`, `:basic`"""

        seed = get(kws, :seed, get(default_kw_args, :seed, 42))

        modular = get(kws, :modular, get(default_kw_args, :modular, :auto))
        @assert modular in (:auto, :classic_modular, :learn_and_apply) """
        Not recognized modular strategy: $modular
        Possible choices for keyword "modular" are:
        `:auto`, `:classic_modular`, `:learn_and_apply`"""

        _composite = get(kws, :_composite, get(default_kw_args, :_composite, 4))

        selection = get(kws, :selection, get(default_kw_args, :selection, :auto))
        @assert selection in (:auto, :normal, :sugar, :be_divided_and_perish)

        check = get(kws, :check, get(default_kw_args, :check, true))
        homogenize = get(kws, :homogenize, get(default_kw_args, :homogenize, :auto))
        @assert homogenize in (:auto, :no, :yes) """
        Not recognized homogenization strategy: $homogenize
        Possible choices for keyword "homogenize" are:
        `:auto`, `:no`, `:yes`"""

        changematrix = get(kws, :changematrix, get(default_kw_args, :changematrix, false))

        _generic = get(kws, :_generic, get(default_kw_args, :_generic, false))

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
            selection,
            modular,
            _composite,
            check,
            homogenize,
            changematrix,
            _generic
        )
    end
end
