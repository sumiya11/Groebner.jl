# Handling keyword arguments in the interface.

#! format: off
# Syntax formatting is disabled for the next couple of lines.
#
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
        :loglevel => _default_logging_level,
        :maxpairs => 0
    ),
    :normalform => Base.ImmutableDict(
        :check    => true,
        :ordering => InputOrdering(),
        :loglevel => _default_logging_level
    ),
    :isgroebner => Base.ImmutableDict(
        :ordering => InputOrdering(),
        :loglevel => _default_logging_level
    ),
    :kbase => Base.ImmutableDict(
        :check    => true,
        :ordering => InputOrdering(),
        :loglevel => _default_logging_level
    )
)
#! format: on

"""
    KeywordsHandler

Stores keyword arguments passed to one of the functions in the interface. On
creation, checks that the arguments are correct.

Sets the global logging level.
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

    function KeywordsHandler(function_key, kws)
        @assert function_key in _supported_kw_arguments
        default_kw_args = _supported_kw_args[function_key]
        @assert all(in(default_kw_args), kw)
        reduced = get(kws, :reduced, get(default_kw_args, :reduced, true))
        ordering = get(kws, :ordering, get(default_kw_args, :ordering, InputOrdering()))
        certify = get(kws, :certify, get(default_kw_args, :certify, false))
        linalg = get(kws, :linalg, get(default_kw_args, :linalg, :prob))
        @assert linalg in (:prob, :det)
        monoms = get(kws, :monoms, get(default_kw_args, :monoms, default_monom_representation()))
        seed = get(kws, :seed, get(default_kw_args, :seed, 42))
        loglevel = get(kws, :loglevel, get(default_kw_args, :loglevel, 0))
        update_logging_level(LogLevel(loglevel))
        # if Threads.threadid() != 1
        #     @log LogLevel(1000) "It looks like Groebner is run in parallel. Note that logging is disabled for all threads except the main thread."
        # end
        maxpairs = get(kws, :linalg, get(default_kw_args, :maxpairs, 0))
        @assert maxpairs >= 0
        new(reduced, ordering, certify, linalg, monoms, seed, loglevel, maxpairs)
    end
end
