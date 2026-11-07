# Prepare once; benchmark workers only activate this resolved environment.
using Pkg
root = normpath(joinpath(@__DIR__, "../../../.."))
environment = joinpath(root, ".local", "environment")
mkpath(environment)
project = joinpath(environment, "Project.toml")
if !isfile(project)
    cp(joinpath(@__DIR__, "../../Project.toml"), project)
    open(project, "a") do io
        println(io, "\n[compat]\nPrettyTables = \"2\"")
    end
end
Pkg.activate(environment)
Pkg.develop(path=root)
Pkg.add(PackageSpec(name="PrettyTables", version="2"))
Pkg.instantiate()
Pkg.precompile()
println("BENCHMARK_JULIA_ENV=$environment")
