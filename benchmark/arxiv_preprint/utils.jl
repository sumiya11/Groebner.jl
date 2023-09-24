# shared utilities for benchmarking
using Printf

const BENCHMARK_RESULTS = "results"

function get_benchmark_dir(id)
    "$BENCHMARK_RESULTS/benchmark_$id"
end

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

timings_filename() = generic_filename("timings")
function timings_filename(id)
    generic_filename("timings", id)
end

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
