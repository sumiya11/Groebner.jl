import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()

using CpuId, Logging, Pkg, Printf
using Statistics

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

include("../utils.jl")

const runtime = Dict()

const PROBLEM_NAME = ARGS[1]
const NUM_RUNS = parse(Int, ARGS[2])
const BENCHMARK_SET = parse(Int, ARGS[3])
const LIB_PATH = ARGS[4]
index1 = startswith(LIB_PATH, "/") ? 2 : 1
index2 = endswith(LIB_PATH, "/") ? length(LIB_PATH) - 1 : length(LIB_PATH)
const LIB_PATH_NORM = LIB_PATH[index1:index2]
const BENCHMARK_DIR = "../../" * get_benchmark_dir("openf4", BENCHMARK_SET)

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

    problempath = (@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/"
    cmd_compille = Cmd(
        `g++ $(problempath)$(PROBLEM_NAME).cpp -I/$LIB_PATH_NORM/include/ -L/$LIB_PATH_NORM/lib/ -lopenf4 -Wl,-rpath=/$LIB_PATH_NORM/lib/ -o $(problempath)$(PROBLEM_NAME)`
    )
    proc = run(cmd_compille, wait=true)
    if proc.exitcode != 0
        @warn "Failed to compile for openf4"
        exit(1)
    end
    for iter in 1:NUM_RUNS
        @info "Computing GB.." iter
        cd("$(problempath)")
        cmd_exec = Cmd(["./$PROBLEM_NAME"])
        timing = @timed proc = run(cmd_exec, wait=true)
        @assert process_exited(proc)
        if proc.exitcode != 0
            @warn "Something probably went wrong in openf4"
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
