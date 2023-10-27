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

"""
    performance_counters_enabled() -> Bool

If performance-tracking macro `@timeit` should be enabled in Groebner. 

When this is `false`, all performance counters in Groebner are disabled and
entail no runtime overhead.
"""
performance_counters_enabled() = false

const _groebner_timer = TimerOutputs.TimerOutput()

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

Note: `@timeit` cannot wrap separate expressions that contain `@label` or
`@goto` (functions are fine).
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
        if $(@__MODULE__).performance_counters_enabled()
            $timed_expr
        else
            $expr
        end
    end
end

function refresh_performance_counters!()
    !performance_counters_enabled() && return nothing
    threadid() != 1 && return nothing
    TimerOutputs.reset_timer!(_groebner_timer)
    nothing
end

function log_performance_counters(statistics)
    (statistics === :no) && return nothing
    if !performance_counters_enabled()
        @log level = 0 """
        Performance counters are not printed since `performance_counters_enabled()` is `false`.
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
