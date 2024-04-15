# Test for performance regressions.
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse, GitHubActions, GitHub, Random, Logging
using Test, TestSetExtensions, InteractiveUtils, PrettyTables
using Base.Threads

const MAX_DEVIATION = 0.1
const IGNORE_SMALL = 1e-3

const dir_master = (@__DIR__) * "/run-on-master"
const dir_nightly = (@__DIR__) * "/run-on-nightly"

function runbench()
    @info "Start benchmarking.."
    @info "Using $(nthreads()) threads"

    # Run benchmarks on master
    @info "Benchmarking Groebner.jl, master, running $dir_master"
    @time run(
        `$(Base.julia_cmd()) --startup-file=no --threads=$(nthreads()) --project=$dir_master $dir_master/run_benchmarks.jl`,
        wait=true
    )

    # Run benchmarks on nightly
    @info "Benchmarking Groebner.jl, nightly, running $dir_nightly"
    @time run(
        `$(Base.julia_cmd()) --startup-file=no --threads=$(nthreads()) --project=$dir_nightly $dir_nightly/run_benchmarks.jl`,
        wait=true
    )
end

# Adapted from https://github.com/MakieOrg/Makie.jl/blob/v0.21.0/metrics/ttfp/run-benchmark.jl.
# License is MIT.
function best_unit(m)
    if m < 1e3
        return 1, "ns"
    elseif m < 1e6
        return 1e3, "μs"
    elseif m < 1e9
        return 1e6, "ms"
    else
        return 1e9, "s"
    end
end

# Compare results
function compare()
    results_master = nothing
    results_nightly = nothing
    try
        results_master_file = open(dir_master * "/results", "r")
        @info "Reading $results_master_file"
        results_master = readlines(results_master_file)
        close(results_master_file)
    catch e
        @warn "Error when reading the file with results"
    end
    try
        results_nightly_file = open(dir_nightly * "/results", "r")
        @info "Reading $results_nightly_file"
        results_nightly = readlines(results_nightly_file)
        close(results_nightly_file)
    catch e
        @warn "Error when reading the file with results"
    end
    @assert !(results_master === nothing) && !(results_nightly === nothing)
    @assert length(results_master) == length(results_nightly)
    @assert !isempty(results_master)
    @assert results_master[1] == "master" && results_nightly[1] == "nightly"
    results_master = results_master[2:end]
    results_nightly = results_nightly[2:end]
    table = Matrix{Any}(undef, length(results_master), 4)
    fail = false
    for (i, (master, nightly)) in enumerate(zip(results_master, results_nightly))
        problem_name_master, times_master = split(master, ":")
        problem_name_nightly, times_nightly = split(nightly, ":")
        @assert problem_name_master == problem_name_nightly
        times_master = map(
            x -> 1e9 * parse(Float64, String(strip(x, ['[', ']', ' ', '\t']))),
            split(times_master, ",")
        )
        times_nightly = map(
            x -> 1e9 * parse(Float64, String(strip(x, ['[', ']', ' ', '\t']))),
            split(times_nightly, ",")
        )
        f, unit = best_unit(maximum(times_master))
        m1 = round(mean(times_master) / f, digits=2)
        d1 = round(std(times_master) / f, digits=2)
        label_master = "$m1 ± $d1 $unit"
        m2 = round(mean(times_nightly) / f, digits=2)
        d2 = round(std(times_nightly) / f, digits=2)
        label_nightly = "$m2 ± $d2 $unit"
        indicator = if mean(times_master) < 1e9 * IGNORE_SMALL
            0, "insignificant"
        elseif (1 + MAX_DEVIATION) * m1 < m2
            fail = true
            2, "**slower**❌"
        elseif m1 > (1 + MAX_DEVIATION) * m2
            1, "**faster**✅"
        else
            0, "insignificant"
        end
        table[i, 1] = problem_name_master
        table[i, 2] = label_master
        table[i, 3] = label_nightly
        table[i, 4] = indicator[2]
    end
    fail, table
end

function post(fail, table)
    comment_header = """
    ## Running times benchmark

    Note, that these numbers may fluctuate on the CI servers, so take them with a grain of salt.

    """
    io = IOBuffer()
    println(io, comment_header)
    if fail
        println(io, "Potential regressions detected❌")
    else
        println(io, "No regressions detected✅")
    end
    table_header = ["Problem", "Master", "This commit", "Result"]
    pretty_table(io, table, header=table_header)
    comment_str = String(take!(io))
    println(comment_str)
end

function main()
    runbench()
    fail, table = compare()
    post(fail, table)
    @test !fail
    versioninfo(verbose=true)
end

main()
