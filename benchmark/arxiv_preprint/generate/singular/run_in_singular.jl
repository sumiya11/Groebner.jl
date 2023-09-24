import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()

using CpuId, Logging, Pkg, Printf
using Statistics

using Singular

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

include("../utils.jl")

const runtime = Dict()

const PROBLEM_NAME = ARGS[1]
const NUM_RUNS = parse(Int, ARGS[2])
const BENCHMARK_SET = parse(Int, ARGS[3])
const BENCHMARK_DIR = "../../" * get_benchmark_dir("singular", BENCHMARK_SET)

@info "" ARGS
@info "" PROBLEM_NAME
@info "" NUM_RUNS
@info "" BENCHMARK_SET
flush(stdout)
flush(stderr)

# Load the system
path = (@__DIR__) * "/$BENCHMARK_DIR/$PROBLEM_NAME/$PROBLEM_NAME.jl"
include(path)

function aa_system_to_singular(system)
    R = AbstractAlgebra.parent(system[1])
    modulo = AbstractAlgebra.characteristic(R)
    n = AbstractAlgebra.nvars(R)
    ground_s = Singular.N_ZpField(modulo)
    R_s, _ = Singular.PolynomialRing(ground_s, ["x$i" for i in 1:n], ordering=:degrevlex)
    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
            ground_s,
            AbstractAlgebra.map_coefficients(c -> ground_s(c.d), f),
            parent=R_s
        ),
        system
    )
    ideal_s = Singular.Ideal(R_s, system_s)
    ideal_s
end

# Compile
singular_system = aa_system_to_singular(system)
Singular.std(system, complete_reduction=true)

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
        singular_system = aa_system_to_singular(system)
        timing = @timed result = Singular.std(system, complete_reduction=true)
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
