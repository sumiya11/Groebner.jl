# shared utilities for benchmarking
using Printf

# Benchmark results directory
const BENCHMARK_RESULTS = "results"
# Directory with data necessary for validating benchmark results
const BENCHMARK_VALIDATE = "validate"
# Table with results is written to file..
const BENCHMARK_TABLE = "benchmark_result"

# Benchmark directory for a specific backend and benchmark set id
function get_benchmark_dir(backend, id)
    "$BENCHMARK_RESULTS/$backend/benchmark_$id"
end

function get_validate_dir(id)
    "$BENCHMARK_VALIDATE/benchmark_$id"
end

# Statistics of interest
const TIME_CATEGORIES = [:total_time]
const DATA_CATEGORIES = []
const ALL_CATEGORIES = union(TIME_CATEGORIES, DATA_CATEGORIES)

const HUMAN_READABLE_CATEGORIES = Dict(:total_time => "Total")
const CATEGORY_FORMAT = Dict()
for cat in ALL_CATEGORIES
    CATEGORY_FORMAT[cat] = (val) -> if val isa Real
        @sprintf("%.2f", val)
    else
        string(val)
    end
end

# Code templates
function julia_pkg_preamble(dir)
    Pkg.Registry.add("General")
    Pkg.activate(dir)
    Pkg.resolve()
    Pkg.instantiate()
    nothing
end

timings_filename() = generic_filename("timings")
function timings_filename(id)
    generic_filename("timings", id)
end

certificate_filename() = generic_filename("certificate")
output_filename() = generic_filename("output")

function result_filename(id)
    generic_filename("result", id)
end

function data_filename(id)
    generic_filename("data", id)
end

function generic_filename(name, id)
    "$(name)_$(id)"
end

function generic_filename(name)
    "$(name)"
end
