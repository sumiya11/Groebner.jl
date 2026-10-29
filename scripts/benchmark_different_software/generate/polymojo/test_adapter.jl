using Test
import AbstractAlgebra, Groebner
include("generate.jl")
include("validate.jl")

function reference_output(system)
    ring = parent(first(system))
    options = isdefined(Base, :pkgversion) && Base.pkgversion(Groebner) >= v"0.10" ? (; tasks=1) : (; threaded=:no)
    result = Groebner.groebner(system; options...)
    join(string.(AbstractAlgebra.gens(ring)), ", ") * "\n" *
        string(AbstractAlgebra.characteristic(AbstractAlgebra.base_ring(ring))) * "\n" *
        join(string.(result), ",\n") * "\n"
end

@testset "Poly Mojo interchange and certificate modes" begin
    for field in (AbstractAlgebra.QQ, AbstractAlgebra.GF(7), AbstractAlgebra.GF(1031), AbstractAlgebra.GF(1073741827))
        _, (x, y) = AbstractAlgebra.polynomial_ring(field, ["x", "y"], internal_ordering=:degrevlex)
        system = [x^2 + y^2 - 1, x + y]
        for mode in ("plain", "certified")
            mktempdir() do dir
                generate_benchmark_file_for_polymojo(system, dir)
                run(Cmd([ENV["POLY_MOJO_DRIVER"], joinpath(dir, "input"), dir,
                         "dummy", "2", mode, "true"]))
                @test compare_full_bases(reference_output(system), read(joinpath(dir, "output"), String))
                @test length(readlines(joinpath(dir, "samples.csv"))) == 3
                @test occursin("verification_time", read(joinpath(dir, "timings"), String)) == (mode == "certified")
            end
        end
    end
    _, (a, b, c) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["a", "b", "c"], internal_ordering=:degrevlex)
    system = [a + AbstractAlgebra.QQ(big(2)^100 + 17, 37)*b, b - c, c - 1]
    for mode in ("plain", "certified")
        mktempdir() do dir
            generate_benchmark_file_for_polymojo(system, dir)
            run(Cmd([ENV["POLY_MOJO_DRIVER"], joinpath(dir, "input"), dir,
                     "large rational", "2", mode, "true"]))
            @test compare_full_bases(reference_output(system), read(joinpath(dir, "output"), String))
            # Truncation must fail instead of producing an apparently valid basis.
            text = read(joinpath(dir, "input"), String)
            write(joinpath(dir, "input"), text[1:end-2])
            proc = run(pipeline(ignorestatus(Cmd([ENV["POLY_MOJO_DRIVER"], joinpath(dir, "input"), dir,
                                                 "bad input", "1", mode, "true"])), stdout=devnull, stderr=devnull))
            @test proc.exitcode != 0
        end
    end
    # Runtime variable counts exercise both packed widths and the generic path.
    for nvars in (9, 17)
        _, variables = AbstractAlgebra.polynomial_ring(AbstractAlgebra.GF(1031),
            ["x$i" for i in 1:nvars], internal_ordering=:degrevlex)
        system = [variables[i] - i for i in 1:nvars]
        for mode in ("plain", "certified")
            mktempdir() do dir
                generate_benchmark_file_for_polymojo(system, dir)
                run(Cmd([ENV["POLY_MOJO_DRIVER"], joinpath(dir, "input"), dir,
                         "many variables", "2", mode, "true"]))
                @test compare_full_bases(reference_output(system), read(joinpath(dir, "output"), String))
            end
        end
    end
end
