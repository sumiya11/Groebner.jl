import Pkg

include("../utils.jl")
julia_pkg_preamble("$(@__DIR__)")

using CpuId, Logging, Pkg, Printf
using Statistics

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

const runtime = Dict()

const PROBLEM_NAME = ARGS[1]
const NUM_RUNS = parse(Int, ARGS[2])
const BENCHMARK_SET = parse(Int, ARGS[3])
const VALIDATE = parse(Bool, ARGS[4])
const BIN_PATH = ARGS[5]
index1 = startswith(BIN_PATH, "/") ? 1 : 1
index2 = endswith(BIN_PATH, "/") ? length(BIN_PATH) - 1 : length(BIN_PATH)
const BIN_PATH_NORM = BIN_PATH[index1:index2]

const BENCHMARK_DIR = "../../" * get_benchmark_dir("maple", BENCHMARK_SET)

@info "" ARGS
@info "" PROBLEM_NAME
@info "" NUM_RUNS
@info "" BENCHMARK_SET
@info "" VALIDATE
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
        problemfile = (@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$(PROBLEM_NAME).mpl"
        cmd = Cmd([BIN_PATH_NORM, problemfile])

        timing = time_ns()
        proc = run(cmd, wait=true)
        timing = (time_ns() - timing) / 1e9
        @assert process_exited(proc)
        if proc.exitcode != 0
            @warn "Something probably went wrong in maple"
            exit(1)
        end
        @debug "Result is" result
        runtime[PROBLEM_NAME][:total_time] = min(runtime[PROBLEM_NAME][:total_time], timing)
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
# dump_timings()
