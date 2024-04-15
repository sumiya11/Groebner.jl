# Test for performance regressions.
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse, GitHubActions, GitHub, Random, Logging
using Test, TestSetExtensions, InteractiveUtils
using Base.Threads

const MAX_ACCEPTABLE_RELATIVE_DEVIATION = 0.1
const IGNORE_SMALL_ABSOLUTE_DEVIATION = 1e-3

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
    @info "Benchmarking Groebner.jl, nightly, running $dir_nightly" dir_nightly
    @time run(
        `$(Base.julia_cmd()) --startup-file=no --threads=$(nthreads()) --project=$dir_nightly $dir_nightly/run_benchmarks.jl`,
        wait=true
    )
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
        problem_name_nightly, times_nightly = split(nightly, ",")
        @assert problem_name_master == problem_name_nightly
        times_master = map(
            x -> parse(Float64, String(strip(x, ['[', ']', ' ', '\t']))),
            split(times_master, ",")
        )
        times_nightly = map(
            x -> parse(Float64, String(strip(x, ['[', ']', ' ', '\t']))),
            split(times_nightly, ",")
        )
        label_master = "$(mean(times_master)) ± $(std(time_master))"
        label_nightly = "$(mean(times_nightly)) ± $(std(times_nightly))"
        indicator =
            if (1 + MAX_ACCEPTABLE_RELATIVE_DEVIATION) * mean(times_master) <
               mean(times_nightly)
                fail = true
                2, "**slower**❌"
            elseif mean(times_master) >
                   (1 + MAX_ACCEPTABLE_RELATIVE_DEVIATION) * mean(times_nightly)
                1, "**faster**✅"
            else
                0, "insignificant"
            end
        table[i, 1] = problem_name
        table[i, 2] = label_master
        table[i, 3] = label_nightly
        table[i, 4] = indicator[2]
        # delta = time_nightly - time_master
        # if abs(delta) / time_master > MAX_ACCEPTABLE_RELATIVE_DEVIATION
        #     if abs(delta) > IGNORE_SMALL_ABSOLUTE_DEVIATION
        #         if delta > 0
        #             printstyled("Regression\n", color=:red)
        #         else
        #             printstyled("Improvement\n", color=:green)
        #         end
        #         @test delta / time_master < MAX_ACCEPTABLE_RELATIVE_DEVIATION
        #     else
        #         printstyled("Insignificant\n")
        #         @test true
        #     end
        # else
        #     printstyled("Insignificant\n")
        #     @test true
        # end
    end
    fail, table
end

# Adapted from https://github.com/MakieOrg/Makie.jl/blob/v0.21.0/metrics/ttfp/run-benchmark.jl.
# License is MIT.
function github_context()
    owner = "sumiya11"
    return (
        owner=owner,
        repo=GitHub.Repo("$(owner)/Groebner.jl"),
        auth=GitHub.authenticate(ENV["GITHUB_TOKEN"])
    )
end

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

function make_or_edit_comment(ctx, pr, benchmarks)
    prev_comments, _ = GitHub.comments(ctx.repo, pr; auth=ctx.auth)
    idx = findfirst(c -> c.user.login == "MakieBot", prev_comments)
    if isnothing(idx)
        comment = update_comment(COMMENT_TEMPLATE, package_name, benchmarks)
        println(comment)
        GitHub.create_comment(ctx.repo, pr; auth=ctx.auth, params=Dict("body" => comment))
    else
        old_comment = prev_comments[idx].body
        comment = update_comment(old_comment, package_name, benchmarks)
        println(comment)
        GitHub.edit_comment(
            ctx.repo,
            prev_comments[idx],
            :pr;
            auth=ctx.auth,
            params=Dict("body" => comment)
        )
    end
end

function post(fail, table)
    header = """
    ## Running times benchmark

    Note, that these numbers may fluctuate on the CI servers, so take them with a grain of salt.
    
    """
    io = IOBuffer()
    println(io, header)
    if fail
        println(io, "Potential regressions detected❌")
    else
        println(io, "No regressions detected✅")
    end
    pretty_table(io, table)
    comment_str = String(take!(io))
    println(comment_str)

    # Post to github
    ctx = github_context()
    pr_to_comment = get(ENV, "PR_NUMBER", nothing)
    if !isnothing(pr_to_comment)
        pr = GitHub.pull_request(ctx.repo, pr_to_comment)
        make_or_edit_comment(ctx, pr, table)
    else
        @info "Not commenting, no PR found"
        println(update_comment(COMMENT_TEMPLATE, table))
    end
end

function main()
    runbench()
    fail, table = compare()
    @test !fail
    post(fail, table)
    versioninfo(verbose=true)
end

main()
