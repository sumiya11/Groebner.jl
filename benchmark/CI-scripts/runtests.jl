# Test for performance regressions.

using Test
using TestSetExtensions
using InteractiveUtils

const MAX_ACCEPTABLE_RELATIVE_DEVIATION = 0.1
const IGNORE_SMALL_ABSOLUTE_DEVIATION = 1e-3

@info "Start benchmarking.."

# Run benchmarks on the latest stable version of Groebner.jl
dir_stable = (@__DIR__) * "/run-on-stable"
@info "Benchmarking Groebner.jl, stable" dir_stable
@time run(`julia --project=$dir_stable $dir_stable/run_benchmarks.jl`, wait=true)

# Run benchmarks on the nightly version of Groebner.jl
dir_nightly = (@__DIR__) * "/run-on-nightly"
@info "Benchmarking Groebner.jl, nightly" dir_nightly
@time run(`julia --project=$dir_nightly $dir_nightly/run_benchmarks.jl`, wait=true)

# Compare results
function compare()
    results_stable = nothing
    results_nightly = nothing
    try
        results_stable_file = open(dir_stable * "/results", "r")
        results_stable = readlines(results_stable_file)
        close(results_stable_file)
    catch e
        @warn "Error when reading the file with results, Groebner.jl stable"
    end
    try
        results_nightly_file = open(dir_nightly * "/results", "r")
        results_nightly = readlines(results_nightly_file)
        close(results_nightly_file)
    catch e
        @warn "Error when reading the file with results, Groebner.jl nightly"
    end
    @assert !(results_stable === nothing) && !(results_nightly === nothing)
    @assert length(results_stable) == length(results_nightly)
    @assert !isempty(results_stable)
    @assert results_stable[1] == "stable" && results_nightly[1] == "nightly"
    results_stable = results_stable[2:end]
    results_nightly = results_nightly[2:end]
    println("Groebner.jl, stable v nightly:")
    for (stable, nightly) in zip(results_stable, results_nightly)
        problem_name..., time_stable = split(stable, ",")
        problem_name_..., time_nightly = split(nightly, ",")
        @assert problem_name == problem_name_
        time_stable = parse(Float64, time_stable)
        time_nightly = parse(Float64, time_nightly)
        print(join(problem_name, ","), ":\t$time_stable v $time_nightly s =>\t")
        delta = time_nightly - time_stable
        if abs(delta) / time_stable > MAX_ACCEPTABLE_RELATIVE_DEVIATION
            if abs(delta) > IGNORE_SMALL_ABSOLUTE_DEVIATION
                if delta > 0
                    printstyled("Regression\n", color=:red)
                else
                    printstyled("Improvement\n", color=:green)
                end
                @test delta / time_stable < MAX_ACCEPTABLE_RELATIVE_DEVIATION
            else
                printstyled("Insignificant\n")
                @test true
            end
        else
            printstyled("Insignificant\n")
            @test true
        end
    end
end

@testset "Compare performance" verbose = true begin
    compare()
    versioninfo(verbose=true)
end
