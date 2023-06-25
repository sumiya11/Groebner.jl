# Handling keyword arguments in the interface.
#
# Functions exported by Groebner have similar keyword arguments with largely
# coinciding default values. For that reason, we centrally handle keyword
# arguments here.

_default_kw_check()    = true
_default_kw_loglevel() = 0

#! format: off
# Syntax formatting is disabled for the next couple of lines.

# Maps a function name to a set of supported keyword arguments with their
# default values.
const _supported_kw_args = Base.ImmutableDict(
    :groebner => Base.ImmutableDict(
        :reduced  => true, 
        :ordering => InputOrdering(),
        :certify  => false,
        :linalg   => :prob,
        :monoms   => Best(),
        :seed     => 42,
        :loglevel => _default_kw_loglevel(),
        :maxpairs => 0
    ),
    :normalform => Base.ImmutableDict(
        :check    => _default_kw_check(),
        :ordering => InputOrdering(),
        :loglevel => _default_kw_loglevel()
    ),
    :isgroebner => Base.ImmutableDict(
        :ordering => InputOrdering(),
        :loglevel => _default_kw_loglevel()
    ),
    :kbase => Base.ImmutableDict(
        :check    => _default_kw_check(),
        :ordering => InputOrdering(),
        :loglevel => _default_kw_loglevel()
    )
)
#! format: on

"""
    Keywords

Stores keyword arguments passed to one of the functions in the interface. On
creation, Checks that the keyword arguments are correct.
"""
struct KeywordsHandler
    reduced::Bool
    ordering::Ord
    certify::Bool
    linalg::Symbol
    monoms::M
    seed::Int
    loglevel::Int
    maxpairs::Int

    function Keywords(function_key, kws)
        @assert function_key in _supported_kw_arguments
        default_kw_args = _supported_kw_args[function_key]
        @assert all(in(default_kw_args), kw)
        reduced = get(kws, :reduced, get(default_kw_args, :reduced, true))
        ordering = get(kws, :ordering, get(default_kw_args, :reduced, InputOrdering()))
        certify = get(kws, :certify, get(default_kw_args, :reduced, false))
        linalg = get(kws, :linalg, get(default_kw_args, :reduced, :prob))
        @assert linalg in (:prob, :det)
        monoms = get(kws, :monoms, get(default_kw_args, :reduced, Best()))
        seed = get(kws, :seed, get(default_kw_args, :reduced, 42))
        loglevel = get(kws, :linalg, get(default_kw_args, :reduced, 0))
        @assert loglevel in (0, 1, 2, 3)
        update_logging_level(LogLevel(loglevel))
        maxpairs = get(kws, :linalg, get(default_kw_args, :reduced, 0))
        @assert maxpairs >= 0
        new(reduced, ordering, certify, linalg, monoms, seed, loglevel, maxpairs)
    end
end
