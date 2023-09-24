using CpuId, Logging, Pkg, Printf
using Statistics

using Groebner

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

include("utils.jl")

const runtime = Dict()

const PROBLEM_NAME = ARGS[1]
const NUM_RUNS = parse(Int, ARGS[2])
const BENCHMARK_SET = parse(Int, ARGS[3])
const BENCHMARK_DIR = get_benchmark_dir(BENCHMARK_SET)

@info "" ARGS
@info "" PROBLEM_NAME
@info "" NUM_RUNS
@info "" BENCHMARK_SET
flush(stdout)
flush(stderr)

# Load the system
path = (@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$PROBLEM_NAME.jl"
include(path)

# Compile
groebner(system)

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
        timing = @timed result = groebner(system)
        @debug "Result is" result
        # for cat in ID_TIME_CATEGORIES
        #     if haskey(StructuralIdentifiability._runtime_logger, cat)
        #         runtime[PROBLEM_NAME][cat] = StructuralIdentifiability._runtime_logger[cat]
        #     end
        # end
        # for cat in ID_runtime_CATEGORIES
        #     if haskey(StructuralIdentifiability._runtime_logger, cat)
        #         runtime[PROBLEM_NAME][cat] =
        #             deepcopy(StructuralIdentifiability._runtime_logger[cat])
        #     end
        # end
        runtime[PROBLEM_NAME][:total_time] =
            min(runtime[PROBLEM_NAME][:total_time], timing.time)
    end
    # for cat in ID_TIME_CATEGORIES
    #     if haskey(runtime[PROBLEM_NAME], cat)
    #         runtime[PROBLEM_NAME][cat] = runtime[PROBLEM_NAME][cat] / NUM_RUNS
    #     end
    # end
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

function dump_results()
    filename = result_filename(GLOBAL_ID)
    open((@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$filename", "w") do io
        if haskey(runtime, PROBLEM_NAME)
            println(io, runtime[PROBLEM_NAME][:return_value])
        end
    end
    filename = data_filename(GLOBAL_ID)
    open((@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$filename", "w") do io
        write(io, "$PROBLEM_NAME\n")
    end
    for cat in ID_DATA_CATEGORIES
        if !haskey(runtime[PROBLEM_NAME], cat)
            continue
        end
        if cat === :something_important
            # make a separate file for it
            filename_cat = generic_filename(cat, GLOBAL_ID)
            open((@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$filename_cat", "w") do io
                # print something
            end
            continue
        end
        # otherwise, print in the data file
        open((@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$filename", "a+") do io
            write(io, "$cat, ")
            write(io, string(runtime[PROBLEM_NAME][cat]))
            write(io, "\n")
        end
    end
end

process_system()
dump_timings()
# dump_results()
