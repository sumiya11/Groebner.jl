# Measuring performance

# Provides the macro @timeit, which can be used to track the performance of code
# blocks and functions. Internally, it uses the TimerOutputs.jl package.
#
# When active, one hit of @timeit has runtime overhead of around 250
# nanoseconds.
# 
# @timeit can be turned off with `performance_counters_enabled`. This eliminates
# all runtime overhead. Thus, it is fine to leave the @timeit statements in the
# source code, since they do not affect performance.

"""
    performance_counters_enabled() -> Bool

If performance-tracking macro `@timeit` should be enabled in Groebner. If
`false`, then all performance counters are disabled, and entail no runtime
overhead.
"""
performance_counters_enabled() = true

timeit_debug_enabled() = false

const _groebner_timer = TimerOutputs.TimerOutput()

"""
    @timeit label expr
    @timeit function foo() ... end

Wraps the expression (or function) in a timed block.

Timed block records some runtime performance information, if performance
counters are enabled in Groebner (see `performance_counters_enabled`).

The function `log_performance_counters` prints the table with statistics into
the current logging IO-stream.

## Example

Assume that `foo` is a function in Groebner.jl that returns the answer to the
question of life, the universe, and everything. Then, we can register it for
benchmarking with

```jldoctest
@timeit function foo()
    42
end
```

We can also benchmark separate expressions:

```jldoctest
@timeit function foo()
    @timeit "initialize data" X = rand(100)
    @timeit "compute stuff" y = map(x -> sin(x) + cos(x), X)
    return y 
end
```
"""
macro timeit(args...)
    esc(quote
        if $(@__MODULE__).performance_counters_enabled()
            Base.@__doc__ TimerOutputs.@timeit($_groebner_timer, $(args...))
        else
            $(args[end])
        end
    end)
end

function refresh_performance_counters!()
    !performance_counters_enabled() && return nothing
    threadid() != 1 && return nothing
    TimerOutputs.reset_timer!(_groebner_timer)
    nothing
end

function log_performance_counters(statistics)
    (statistics === :no) && return nothing
    !performance_counters_enabled() && return nothing
    threadid() != 1 && return nothing
    iostream = _default_logger[].stream
    TimerOutputs.print_timer(
        iostream,
        _groebner_timer,
        allocations=true,
        sortby=:time,
        linechars=:ascii,
        compact=false,
        title="Groebner.jl"
    )
end
