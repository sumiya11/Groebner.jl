using Test
include("../../one_script_to_run_them_all.jl")

@testset "Optional backend and certificate reporting" begin
    @test "polymojo" in _available_backends
    @test "polymojocertified" in _available_backends
    @test length(comparison_columns(["groebner", "singular"])) == 2
    @test length(comparison_columns(["groebner", "polymojocertified"])) == 3
    small_prime_suite = get_benchmark_suite(2)
    @test all("chandra $i" in first.(small_prime_suite.systems) for i in 11:14)
    @test all(AbstractAlgebra.characteristic(base_ring(parent(first(system)))) == 1031
              for (_, system) in small_prime_suite.systems)
    args = Dict{String, Any}("backend" => "polymojo", "suite" => 0,
        "bpolymojo" => ENV["POLY_MOJO_DRIVER"], "nruns" => 2,
        "nworkers" => 1, "timeout" => 60, "validate" => "no")
    @test isnothing(check_args(args))
    invalid = merge(args, Dict("bpolymojo" => ""))
    @test_throws ErrorException check_args(invalid)
    # Other backends need no Mojo setting at all.
    other = Dict{String, Any}("backend" => "groebner", "suite" => 0)
    @test isnothing(check_args(other))

    command = get_command_to_run_benchmark("polymojocertified", "dummy 1", 2, 0, true, args)
    @test command.exec[1] == abspath(ENV["POLY_MOJO_DRIVER"])
    @test command.exec[4:end] == ["dummy 1", "2", "certified", "true"]

    # These outputs come from the upstream suite-0 CLI smoke runs.
    names = ["dummy 1", "dummy 2"]
    runtimes = Dict()
    for backend in ("polymojo", "polymojocertified")
        current = merge(args, Dict("backend" => backend))
        runtimes[backend] = collect_timings(current, names)
        @test all(haskey(runtimes[backend][name], :total_time) for name in names)
        report = read(joinpath(@__DIR__, "../../results", backend, "benchmark_result_0.md"), String)
        @test occursin("Verify certificate, s", report) == (backend == "polymojocertified")
    end
    combined = joinpath(@__DIR__, "../../results/benchmark_result_0.md")
    original = isfile(combined) ? read(combined) : nothing
    try
        collect_all_timings(args, runtimes, names)
        report = read(combined, String)
        @test occursin("polymojocertified (verify)", report)
        expected = CATEGORY_FORMAT[:verification_time](runtimes["polymojocertified"]["dummy 1"][:verification_time])
        @test occursin(expected, report)
        # A timed-out/missing case must display missing cells, not throw.
        delete!(runtimes["polymojocertified"], "dummy 2")
        collect_all_timings(args, runtimes, names)
        row = only(filter(line -> startswith(line, "|dummy 2|"), readlines(combined)))
        @test occursin("| - | - |", row)
        collect_all_timings(args, Dict("polymojo" => runtimes["polymojo"]), names)
        @test !occursin("verify", read(combined, String))
    finally
        if original === nothing
            rm(combined, force=true)
        else
            write(combined, original)
        end
    end
end
