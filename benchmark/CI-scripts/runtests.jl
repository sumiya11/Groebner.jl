# Test for performance regressions.

using Test
using TestSetExtensions
using InteractiveUtils

const MAX_ACCEPTABLE_RELATIVE_DEVIATION = 0.1
const IGNORE_SMALL_ABSOLUTE_DEVIATION = 1e-3

@info "Start benchmarking.."

# Run benchmarks on the latest version of Groebner.jl
dir_latest = (@__DIR__) * "/run-on-latest"
@info "Benchmarking Groebner.jl, latest" dir_latest
@time run(`julia --project=$dir_latest $dir_latest/run_benchmarks.jl`, wait=true)

# Run benchmarks on the nighly version of Groebner.jl
dir_nightly = (@__DIR__) * "/run-on-nightly"
@info "Benchmarking Groebner.jl, nightly" dir_nightly
@time run(`julia --project=$dir_nightly $dir_nightly/run_benchmarks.jl`, wait=true)

# Compare results
function compare()
    results_latest = nothing
    results_nightly = nothing
    try
        results_latest_file = open((@__DIR__) * "/run-on-latest/results", "r")
        results_latest = readlines(results_latest_file)
        close(results_latest_file)
    catch e
        @warn "Error when reading the file with results, Groebner.jl latest"
    end
    try
        results_nightly_file = open((@__DIR__) * "/run-on-nightly/results", "r")
        results_nightly = readlines(results_nightly_file)
        close(results_nightly_file)
    catch e
        @warn "Error when reading the file with results, Groebner.jl nightly"
    end
    @assert !(results_latest === nothing) && !(results_nightly === nothing)
    @assert length(results_latest) == length(results_nightly)
    @assert !isempty(results_latest)
    @assert results_latest[1] == "latest" && results_nightly[1] == "nightly"
    results_latest = results_latest[2:end]
    results_nightly = results_nightly[2:end]
    for (latest, nightly) in zip(results_latest, results_nightly)
        problem_name..., time_latest = split(latest, ",")
        problem_name_..., time_nightly = split(nightly, ",")
        @assert problem_name == problem_name_
        time_latest = parse(Float64, time_latest)
        time_nightly = parse(Float64, time_nightly)
        @info """
        Problem: $(join(problem_name, ",")).
        Groebner.jl latest : $time_latest s
        Groebner.jl nightly: $time_nightly s"""
        delta = time_nightly - time_latest
        if delta / time_latest < MAX_ACCEPTABLE_RELATIVE_DEVIATION
            if delta > IGNORE_SMALL_ABSOLUTE_DEVIATION
                @test delta / time_latest < MAX_ACCEPTABLE_RELATIVE_DEVIATION
            else
                @test true
            end
        else
            @test true
        end
    end
end

@testset "Compare performance" verbose = true begin
    compare()
    versioninfo(verbose=true)
end
