import Pkg

include("../utils.jl")
julia_pkg_preamble("$(@__DIR__)")

using CpuId, Logging, Pkg, Printf
using Statistics

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

include("../utils.jl")

const runtime = Dict()

const PROBLEM_NAME = ARGS[1]
const NUM_RUNS = parse(Int, ARGS[2])
const BENCHMARK_SET = parse(Int, ARGS[3])
const BENCHMARK_DIR = "../../" * get_benchmark_dir("msolve", BENCHMARK_SET)

@info "" ARGS
@info "" PROBLEM_NAME
@info "" NUM_RUNS
@info "" BENCHMARK_SET
@info "" "$(@__DIR__)"
flush(stdout)
flush(stderr)

function process_system()
    @info "Processing $PROBLEM_NAME"
    @info """
    Averaging over $NUM_RUNS runs."""

    runtime[PROBLEM_NAME] = Dict{Any, Any}()
    for cat in ALL_CATEGORIES
        runtime[PROBLEM_NAME][cat] = Inf
    end
    for iter in 1:NUM_RUNS
        @info "Computing GB.." iter
        problemfile = (@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$(PROBLEM_NAME).in"
        cmd = Cmd(`msolve -g 2 -l 44 -c 0 -f $problemfile -o /dev/null`)
        timing = @timed proc = run(cmd, wait=true)
        @assert process_exited(proc)
        if proc.exitcode != 0
            @warn "Something probably went wrong in msolve"
            exit(1)
        end
        @debug "Result is" result
        runtime[PROBLEM_NAME][:total_time] =
            min(runtime[PROBLEM_NAME][:total_time], timing.time)
    end
end

function dump_timings()
    timings = ""
    timings *= "$PROBLEM_NAME\n"
    for (key, model_runtime) in runtime
        for c in [:total_time]
            timings *= "$c, "
            timings *= string(model_runtime[c]) * "\n"
        end
    end
    filename = timings_filename()
    open((@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$filename", "w") do io
        write(io, timings)
    end
end

process_system()
dump_timings()
