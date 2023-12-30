# This file is a part of Groebner.jl. License is GNU GPL v2.

# Collecting and displaying statistics about computation.
# This file is WIP!

# Provides the macro `@stat`, which can be used to record a useful value
# (statistic) in runtime and print it after the computation finishes.
#

"""
    statistics_enabled() -> Bool

If statistics-collecting macro `@stat` should be enabled in Groebner. 

When this is `false`, all occurrences of `@stat` in Groebner are disabled and
entail no runtime overhead.
"""
statistics_enabled() = false

###
# GroebnerStatistics and some utilities

_f4_by_iteration_keys = (
    :critical_pairs_deg,
    :critical_pairs_count,
    :new_basis_elements,
    :matrix_block_sizes,
    :matrix_nnz
)

_human_readable_keys = Dict(
    :critical_pairs_deg   => "pair deg.",
    :critical_pairs_count => "# pairs",
    :matrix_block_sizes   => "matrix size",
    :matrix_nnz           => "matrix nnz",
    :new_basis_elements   => "# new pivots"
)

# _f4_numerical_stats = (:f4_iterations,)

mutable struct GroebnerStatistics
    stopwatch_start::Int

    # 
    general_stats::Any

    # Statistics from F4
    F4_iterations::Int
    F4_stats::Any

    # Statistics from multi-modular computation
    multimodular_stats::Any

    other_stats::Any

    function GroebnerStatistics()
        new(0, Dict(), 0, Dict(), Dict(), Dict())
    end
end

function update_statistic!(
    stat::GroebnerStatistics,
    file,
    line,
    args::Any,
    key::Symbol,
    value::Any
)
    if key in _f4_by_iteration_keys
        if !haskey(stat.F4_stats, key)
            stat.F4_stats[key] = Vector{typeof(value)}()
        end
        push!(stat.F4_stats[key], value)
    elseif key === :f4_iterations
        stat.F4_iterations = value
    end

    nothing
end

function refresh_statistics!(stat::GroebnerStatistics)
    stat.stopwatch_start = time_ns()
    empty!(stat.F4_stats)
    nothing
end

function print_f4_statistics(io, stat::GroebnerStatistics)
    F4_iterations = stat.F4_iterations
    if iszero(F4_iterations)
        println(io, "No iterations of F4 were performed!")
        return nothing
    end

    ncols = length(_f4_by_iteration_keys) + 1
    header = Vector{String}(undef, ncols)
    table = Matrix{Any}(undef, F4_iterations, ncols)
    header[1] = ""
    table[:, 1] = collect(1:F4_iterations)
    for (i, key) in enumerate(_f4_by_iteration_keys)
        header[i + 1] = _human_readable_keys[key]
        table[:, i + 1] .= stat.F4_stats[key][1:F4_iterations]
    end
    PrettyTables.pretty_table(io, table, header=header, tf=PrettyTables.tf_compact)

    nothing
end

function print_statistics(io, stat::GroebnerStatistics)
    println(io, "Groebner statistics")

    print_f4_statistics(io, stat)

    nothing
end

const _groebner_statistics = GroebnerStatistics()

function update_statistic(file, line, args, key, value)
    update_statistic!(_groebner_statistics, file, line, args, key, value)
end

function refresh_statistics()
    !statistics_enabled() && return nothing
    threadid() != 1 && return nothing
    refresh_statistics!(_groebner_statistics)
    nothing
end

###
# The @stat macro

function process_stat_args(exprs)
    kwargs = Any[]
    for ex in exprs
        if ex isa Expr && ex.head === :(=)
            k, v = ex.args
            if !(k isa Symbol)
                k = Symbol(k)
            end
            push!(kwargs, (k, v))
        else
            push!(kwargs, (Symbol(ex), ex))
        end
    end
    [], kwargs
end

macro stat(exprs...)
    file, line = String(__source__.file), Int(__source__.line)
    args, kwargs = process_stat_args(exprs)
    body = :()
    for (key, value) in kwargs
        key = Meta.quot(key)
        body = quote
            $body
            update_statistic($file, $line, $args, $key, $value)
        end
    end
    esc(quote
        if $(@__MODULE__).statistics_enabled()
            $body
        else
            nothing
        end
    end)
end

function print_statistics(statistics)
    (statistics in (:no, :timings)) && return nothing
    if statistics in (:stats, :all) && !statistics_enabled()
        @log level = 1_000 """
        Statistics were not collected since `statistics_enabled()` is `false`.
        Consider setting `Groebner.statistics_enabled()` to `true` and trying again.
        """
        return nothing
    end
    threadid() != 1 && return nothing

    print_statistics(stdout, _groebner_statistics)
end
