function generate_benchmark_file_for_polymojo(system, dir)
    ring = parent(first(system))
    field = AbstractAlgebra.base_ring(ring)
    characteristic = AbstractAlgebra.characteristic(field)
    open(joinpath(dir, "input"), "w") do io
        println(io, "POLYMOJO GROEBNER 1")
        println(io, characteristic)
        println(io, AbstractAlgebra.nvars(ring))
        foreach(v -> println(io, v), AbstractAlgebra.gens(ring))
        println(io, length(system))
        for polynomial in system
            println(io, length(polynomial))
            for (coefficient, exponents) in zip(AbstractAlgebra.coefficients(polynomial), AbstractAlgebra.exponent_vectors(polynomial))
                println(io, replace(string(coefficient), "//" => "/"))
                foreach(e -> println(io, e), exponents)
            end
        end
    end
end
