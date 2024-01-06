# This file is a part of Groebner.jl. License is GNU GPL v2.

# Measuring performance

# Provides the macro @timeit, which can be used to track the performance of code
# blocks and functions. Internally, it uses the TimerOutputs.jl package.
#
# When active, one hit of @timeit has runtime overhead of around 250
# nanoseconds.
# 
# @timeit can be turned off by setting `performance_counters_enabled` to
# `false`. This will compile out `@timeit` statements (and eliminate all runtime
# overhead). Thus, it is fine to leave the @timeit statements in the source
# code, since they do not affect performance when inactive.

const _groebner_timer = TimerOutputs.TimerOutput()
const _groebner_timer_lock = Ref{ReentrantLock}(ReentrantLock())

@noinline __throw_timeit_error() = throw(ArgumentError("""
    Invalid usage of macro @timeit in Groebner.jl. Use it as:
    @timeit label expr
    @timeit function foo() ... end"""))

"""
    @timeit label expr
    @timeit function foo() ... end

Wraps the expression (or function) in a timed block.

Timed block records some runtime performance information, if performance
counters are enabled in Groebner (see `performance_counters_enabled`).

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
    @timeit "initialize data" xs = rand(100)
    @timeit "compute stuff" y = map(x -> sin(x) + cos(x), xs)
    return y 
end
```

NOTE: `@timeit` cannot wrap code blocks that contain `@label` or `@goto`.
NOTE: `@timeit` is broken in combination with `Base.Threads.@threads`.
"""
macro timeit(args...)
    if isempty(args) || length(args) > 2
        __throw_timeit_error()
    end
    _timeit(__module__, args...)
end

function _timeit(m, expr)
    _timeit(m, nothing, expr)
end

function _timeit(m, label, expr)
    timed_expr = if TimerOutputs.is_func_def(expr)
        _timeit_func(m, label, expr)
    else
        _timeit_expr(m, label, expr)
    end
    esc(timed_expr)
end

function _timeit_func(m, label, expr)
    expr = macroexpand(m, expr)
    def = splitdef(expr)
    label === nothing && (label = string(def[:name]))
    @debug "Groebner.@timeit now tracks function $label"
    def[:body] = _timeit_expr(m, label, def[:body])
    combinedef(def)
end

function _timeit_expr(m, label, expr)
    # NOTE: some care should be taken to handle the case when expr contains
    # statements such as @label and @goto 
    timed_expr = TimerOutputs._timer_expr(m, false, _groebner_timer, label, expr)
    quote
        if $m.performance_counters_enabled()
            $timed_expr
        else
            $expr
        end
    end
end

function performance_counters_refresh()
    # Do nothing if counters are disabled
    !performance_counters_enabled() && return nothing
    # Do nothing if run from a worker thread.
    # NOTE: this does not always do what is intended. It is still correct, since
    # we lock the timer anyway.
    threadid() != 1 && return nothing

    lock(_groebner_timer_lock[])
    try
        TimerOutputs.reset_timer!(_groebner_timer)
    finally
        unlock(_groebner_timer_lock[])
    end

    nothing
end

function performance_counters_print(statistics)
    (statistics in (:no, :stats)) && return nothing
    if statistics in (:timings, :all) && !performance_counters_enabled()
        @log level = 1_000 """
        Timings were not collected since `performance_counters_enabled()` is `false`.
        Consider setting `Groebner.performance_counters_enabled()` to `true` and trying again.
        """
        return nothing
    end
    threadid() != 1 && return nothing

    iostream = stdout
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
