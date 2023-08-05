# Measuring performance

# Provides the macro @timed_block, which can be used to track the performance of
# a code block / function

"""
    rdtsc() -> Int

Returns the Read Time Stamp Counter.
"""
rdtsc() = ccall("llvm.x86.rdtsc", llvmcall, Int, ())

"""
    performance_counters_enabled() -> Bool

Specifies if performance should be tracked in Groebner. If `false`, then all
performance counters are disabled, and entail no runtime overhead.

See also `@timed_block`.
"""
performance_counters_enabled() = true

"""
    @timed_block begin ... end
    @timed_block function ... end

Wraps the expression in a timed block.

Timed block tracks some runtime performance information, if performance counters
are enabled in Groebner (see `performance_counters_enabled`).

## Example

```jldoctest
@timed_block function dot_product(x, y)
    sum(x .* y)
end

for i in 1:100
    dot_product([1, 2, 3], [4, 5, 6])
end

@log_performance_counters
# prints the following (or something similar)

                            hits    cycles/hit
REPL[1]:1/dot_product       100     91.8
```

*Note that this macro measures performance correctly only for the code that is
defined in the Groebner module. In that sense, the example may be misleading.*
"""
macro timed_block(expr)
    @assert expr.head === :block ||
            expr.head === :function ||
            (expr.head === :(=) && expr.args[1] === :call)
    file, line = String(__source__.file), Int(__source__.line)
    if expr.head === :block
        # if begin ... end
        id = (file, line)
        esc(_timed_block(file, line, id, expr))
    else
        # if function
        fdef = splitdef(expr)
        id = fdef[:name]
        fdef[:body] = _timed_block(file, line, id, fdef[:body])
        esc(combinedef(fdef))
    end
end

"""
    @log_performance_counters

Prints performance counters to the current logging stream.
"""
macro log_performance_counters()
    esc(quote
        if $(@__MODULE__).performance_counters_enabled()
            _log_performance_counters()
        else
            nothing
        end
    end)
end

mutable struct PerfCounterRecord
    id::Any
    file::String
    line::Int
    hits::Int
    cycles::Int
    PerfCounterRecord() = new(nothing, "", 0, 0, 0)
end

const _perf_counters = Vector{PerfCounterRecord}()
const _perf_id_to_index = Dict{Any, Int}()

function _timed_block(file, line, id, expr)
    push!(_perf_counters, PerfCounterRecord())
    _perf_id_to_index[id] = length(_perf_counters)
    index = length(_perf_id_to_index)
    _perf_counters[index].id = id
    _perf_counters[index].file = file
    _perf_counters[index].line = line
    beginstamp = gensym()
    counter = gensym()
    result = gensym()
    quote
        if $(@__MODULE__).performance_counters_enabled()
            $beginstamp = rdtsc()
            $result = $expr
            $counter = rdtsc() - $beginstamp
            accumulate_counter($id, $file, $line, $index, $counter)
            $result
        else
            $expr
        end
    end
end

function accumulate_counter(id, file, line, index, cycles)
    @inbounds _perf_counters[index].hits += 1
    @inbounds _perf_counters[index].cycles += cycles
    nothing
end

function _log_performance_counters()
    io = _default_logger[].stream
    println(io, "\t\t\thits\tcycles/hit")
    for i in 1:length(_perf_counters)
        record = _perf_counters[i]
        println(
            io,
            "$(last(split(record.file, "/"))):$(record.line)/$(record.id)\t\t$(record.hits)\t$(record.cycles/record.hits)"
        )
    end
end
