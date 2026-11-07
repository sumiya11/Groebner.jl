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
const TIME_CATEGORIES = [
    :total_time,
    :total_time_F4,
    :total_time_learn,
    :total_time_apply,
    :total_time_apply_4x,
    :total_time_apply_8x,
    :verification_time
]
const DATA_CATEGORIES = []
const ALL_CATEGORIES = union(TIME_CATEGORIES, DATA_CATEGORIES)

const HUMAN_READABLE_CATEGORIES = Dict(
    :total_time => "Total, s",
    :total_time_learn => "Learn, s",
    :total_time_apply => "Apply, s",
    :total_time_apply_4x => "Apply 4x, s",
    :total_time_apply_8x => "Apply 8x, s",
    :total_time_F4 => "F4, s",
    :verification_time => "Verify certificate, s"
)
const CATEGORY_FORMAT = Dict()
for cat in ALL_CATEGORIES
    CATEGORY_FORMAT[cat] = (val) -> if val isa Real
        @sprintf("%.2f", val)
    else
        string(val)
    end
end

# Verification is supplementary, not part of the comparison's computation time.
CATEGORY_FORMAT[:verification_time] = val -> val isa Real ? @sprintf("%.6g", val) : string(val)

function comparison_columns(backends)
    columns = [(backend=backend, target=:total_time, label=backend) for backend in backends]
    if "polymojocertified" in backends
        push!(columns, (backend="polymojocertified", target=:verification_time,
                        label="polymojocertified (verify)"))
    end
    columns
end

# Code templates
function julia_pkg_preamble(dir)
    # A pre-instantiated environment can be shared by the controller and workers.
    if haskey(ENV, "BENCHMARK_JULIA_ENV")
        Pkg.activate(ENV["BENCHMARK_JULIA_ENV"])
        return nothing
    end
    @info "Activating environment at $dir"
    Pkg.activate(dir)
    Pkg.Registry.add("General")
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
