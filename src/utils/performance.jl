# Measuring performance

# Provides the macro @timed_block, which can be used to track the performance of
# a code block / function

"""
    performance_counters_enabled() -> Bool

Specifies if performance should be tracked in Groebner. If `false`, then all
performance counters are disabled, and entail no runtime overhead.

See also `@timed_block`.
"""
performance_counters_enabled() = true

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
                log_performance_counters()
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
const _perf_stopwatch_start = Ref{UInt}(UInt(0))
const _perf_id_to_index = Dict{Any, Int}()

function refresh_performance_counters!()
    !performance_counters_enabled() && return nothing
    threadid() != 1 && return nothing
    for record in _perf_counters
        refresh_counter!(record)
    end
    _perf_stopwatch_start[] = time_ns()
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
            # $beginstamp = rdtsc()
            $beginstamp = time_ns()
            $result = $expr
            # $counter = rdtsc() - $beginstamp
            $counter = time_ns() - $beginstamp
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

function log_performance_counters(statistics)
    if statistics === :no
        return nothing
    end
    _log_performance_counters()
end

const _perf_table_columns = ["Location", "Hits", "Time / Hit", "Time", "% of total"]

function _log_performance_counters()
    # columns = ["location", "hits", "cycles / hit", "cycles, total"]
    rows = map(counter -> counter.id, _perf_counters)

    m, n = length(rows), length(_perf_table_columns)
    table = Array{Any, 2}(undef, m, n)

    total_cycles_recorded = sum(counter -> counter.cycles, _perf_counters)
    total_cycles = time_ns() - _perf_stopwatch_start[]
    for (i, counter) in enumerate(_perf_counters)
        file = counter.file
        line = counter.line
        hits = counter.hits
        cycles = counter.cycles
        cycles_per_hit = hits == 0 ? 0 : div(cycles, hits)
        percent = cycles / total_cycles
        table[i, 1] = splitpath(file)[end] * ":" * string(line)
        table[i, 2] = hits
        # table[i, 3] = pretty_number(String, cycles_per_hit)
        table[i, 3] = prettytime(cycles_per_hit)
        table[i, 4] = prettytime(cycles)
        table[i, 5] = prettypercent(percent)
    end

    row_permutation = collect(1:m)
    sort!(row_permutation, by=idx -> parse(Float64, table[idx, end][1:(end - 1)]), rev=true)
    table = table[row_permutation, :]
    rows = rows[row_permutation]

    PrettyTables.pretty_table(
        table,
        title="Performance counters for Groebner.jl",
        # label="Does not account for overlap",
        row_labels=rows,
        header=_perf_table_columns,
        limit_printing=false,
        display_size=(-1, -1)
    )
    # println(
    #     "Recorded $(prettytime(total_cycles_recorded)) out of total $(prettytime(total_cycles)) ($(prettypercent(total_cycles_recorded/total_cycles)))."
    # )
    println("Timings do not account for overlap.")
end
