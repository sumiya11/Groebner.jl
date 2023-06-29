# Handling keyword arguments in the interface.

#! format: off
# Syntax formatting is disabled for the next several lines.

# Maps a function name to a set of supported keyword arguments with their
# corresponding default values. default values.
const _supported_kw_args = Dict(
    :groebner => Dict{Symbol, Any}(
        :reduced  => true, 
        :ordering => nothing,
        :certify  => false,
        :linalg   => :randomized,
        :monoms   => nothing,
        :seed     => 42,
        :loglevel => _default_logging_level,
        :maxpairs => 0
    ),
    :normalform => Dict{Symbol, Any}(
        :check    => true,
        :ordering => nothing,
        :loglevel => _default_logging_level
    ),
    :isgroebner => Dict{Symbol, Any}(
        :ordering => nothing,
        :loglevel => _default_logging_level
    ),
    :kbase => Dict{Symbol, Any}(
        :check    => true,
        :ordering => nothing,
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
struct KeywordsHandler{Ord, Monoms}
    reduced::Bool
    ordering::Union{Ord, Nothing}
    certify::Bool
    linalg::Symbol
    monoms::Union{Monoms, Nothing}
    seed::Int
    loglevel::Int
    maxpairs::Int

    function KeywordsHandler(function_key, kws)
        @assert haskey(_supported_kw_args, function_key)
        default_kw_args = _supported_kw_args[function_key]
        @assert all(in(default_kw_args), kws)
        reduced = get(kws, :reduced, get(default_kw_args, :reduced, true))
        ordering = get(kws, :ordering, get(default_kw_args, :ordering, nothing))
        certify = get(kws, :certify, get(default_kw_args, :certify, false))
        linalg = get(kws, :linalg, get(default_kw_args, :linalg, :randomized))
        @assert linalg in (:randomized, :deterministic)
        monoms = get(kws, :monoms, get(default_kw_args, :monoms, nothing))
        seed = get(kws, :seed, get(default_kw_args, :seed, 42))
        loglevel = get(kws, :loglevel, get(default_kw_args, :loglevel, 0))
        update_logging_level(LogLevel(loglevel))
        # if Threads.threadid() != 1
        #     @log LogLevel(1000) "It looks like Groebner is run in parallel. Note that logging is disabled for all threads except the main thread."
        # end
        maxpairs = get(kws, :linalg, get(default_kw_args, :maxpairs, 0))
        @assert maxpairs >= 0
        @log level = 3 """
          Using keywords: 
          reduced=$reduced, 
          ordering=$ordering, 
          certify=$certify, 
          linalg=$linalg, 
          monoms=$monoms, 
          seed=$seed, 
          loglevel=$loglevel, 
          maxpairs=$maxpairs"""
        new{typeof(ordering), typeof(monoms)}(
            reduced,
            ordering,
            certify,
            linalg,
            monoms,
            seed,
            loglevel,
            maxpairs
        )
    end
end
