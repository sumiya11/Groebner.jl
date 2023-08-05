# Measuring performance

# Provides the macro @timed_block, which can be used to track the performance of
# a code block / function

"""
    performance_counters_enabled() -> Bool

Specifies if performance should be tracked in Groebner. If `false`, then all
performance counters are disabled, and entail no runtime overhead.

See also `@timed_block`.
"""
performance_counters_enabled() = false

"""
    rdtsc() -> Int

Returns the Read Time Stamp Counter.
"""
rdtsc() = ccall("llvm.x86.rdtsc", llvmcall, Int, ())

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
        type = :block
        esc(_timed_block(file, line, id, type, expr))
    else
        # if function
        fdef = splitdef(expr)
        id = Meta.quot(fdef[:name])
        type = :function
        fdef[:body] = _timed_block(file, line, id, type, fdef[:body])
        combinedef(fdef)
    end
end

"""
    @log_performance_counters

Logs performance counters to the current logging stream.
"""
macro log_performance_counters()
    esc(
        quote
            if $(@__MODULE__).logging_enabled() &&
               $(@__MODULE__).performance_counters_enabled()
                _log_performance_counters()
            else
                nothing
            end
        end
    )
end

mutable struct PerfCounterRecord
    id::Any
    type::Symbol
    file::String
    line::Int
    hits::Int
    cycles::Int
    PerfCounterRecord() = new(nothing, :none, "", 0, 0, 0)
end

function refresh_counter!(perf::PerfCounterRecord)
    perf.hits = 0
    perf.cycles = 0
end

const _perf_counters = Vector{PerfCounterRecord}()
const _perf_id_to_index = Dict{Any, Int}()

function refresh_performance_counters!()
    for record in _perf_counters
        refresh_counter!(record)
    end
    nothing
end

function _timed_block(file, line, id, type, expr)
    push!(_perf_counters, PerfCounterRecord())
    _perf_id_to_index[id] = length(_perf_counters)
    index = length(_perf_id_to_index)
    _perf_counters[index].id = id
    _perf_counters[index].type = type
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
            accumulate_counter($index, $counter)
            $result
        else
            $expr
        end
    end
end

function accumulate_counter(index, cycles)
    @inbounds _perf_counters[index].hits += 1
    @inbounds _perf_counters[index].cycles += cycles
    nothing
end

function _log_performance_counters()
    msg = ""
    tracked_names = Vector{String}(undef, length(_perf_counters))
    for i in 1:length(_perf_counters)
        record = _perf_counters[i]
        loc = if record.type === :function
            "$(record.id)"
        else
            "$(last(split(record.file, "/"))):($(record.line))"
        end
        tracked_names[i] = loc
    end
    table_left_col_length = maximum(length, tracked_names)
    tab_length = 8  # hmm?
    ntabs = div(table_left_col_length - 1, tab_length) + 2
    header = "Performance counters (does not account for overlap)\n"
    table_header = "\t"^ntabs * "hits\tcycles/hit\n"
    msg *= header
    msg *= table_header
    for i in 1:length(_perf_counters)
        record = _perf_counters[i]
        name = tracked_names[i]
        ntabs_i = 0
        hits = record.hits
        cycles = record.cycles
        mean_cycles = if !iszero(hits)
            round(Int, cycles / hits)
        else
            0
        end
        msg *= name * " "^(tab_length * ntabs - length(name))
        msg *= "\t"^ntabs_i * "$hits\t$mean_cycles\n"
    end
    @log level = -1 msg
    nothing
end
