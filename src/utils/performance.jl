# Helper functions for performance 
# TODO: description,
# TODO: this is not thread-safe

"""
    enable_performance_counters(flag)

Enable debugging globally.
"""
function enable_performance_counters(flag::Bool)
    if !getfield(@__MODULE__, :performance_counters_enabled)()
        Core.eval(@__MODULE__, :(performance_counters_enabled() = $flag))
    end
end

"""
    CounterRecord

"""
mutable struct CounterRecord
    file::String
    line::Int
    hits::Int
    cycles::Int
    PerformanceRecord() = new("", 0, 0, 0)
end

# NOTE: be careful about moving these around!
const __counter_value = Ref{Int}(0)
const __counters = [CounterRecord() for _ in 1:((@__COUNTER__) - 1)]

"""
    @__COUNTER__ -> Int

Expand to an integer that starts at 1. 
The value is incremented by 1 every time it is used in a source file.
"""
macro __COUNTER__()
    __counter_value[] += 1
end

function accumulate_counter(file, line, index, cycles)
    __counters[index].file = file
    __counters[index].line = line
    __counters[index].hits += 1
    __counters[index].cycles += cycles
    if cycles > 5000
        _display_performance_counters()
    end
    nothing
end

"""
    rdtsc() -> Int

Returns rdtsc.

**Note:** the overhead is platform dependent, 
with the average of 10-40 cycles and a high variance.
"""
rdtsc() = ccall("llvm.x86.rdtsc", llvmcall, Int, ())

"""
    @timed_block

The runtime overhead is very small: roughly 20-30 cycles. 
"""
macro timed_block(expr)
    @assert expr.head === :block ||
            expr.head === :function ||
            (expr.head === :(=) && expr.args[1] === :call)
    file, line = String(__source__.file), Int(__source__.line)
    if expr.head === :block
        _timed_block(file, line, expr)
    else # if function
        fdef = splitdef(expr)
        fdef[:body] = _timed_block(file, line, fdef[:body])
        esc(combinedef(fdef))
    end
end

function _timed_block(file, line, expr)
    beginstamp = gensym()
    endstamp = gensym()
    esc(:(
        if $(@__MODULE__).enable_performance_counters()
            $beginstamp = rdtsc()
            $expr
            $endstamp = rdtsc() - $beginstamp
            record_counter($file, $line, @__COUNTER__, $endstamp)
        else
            $expr
        end
    ))
end

macro display_performance_counters()
    esc(:(
        if $(@__MODULE__).enable_performance_counters()
            _display_performance_counters()
        else
        end
    ))
end
function _display_performance_counters()
    println("\t\t\thits\tcycles/hit")
    for i in 1:length(__counters)
        record = __counters[i]
        println(
            "$(last(split(record.file, "/"))):$(record.line)\t\t$(record.hits)\t$(record.cycles/record.hits)"
        )
    end
end
