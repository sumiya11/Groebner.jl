using AbstractAlgebra, Groebner, BenchmarkTools, PrettyTables, Printf, TimerOutputs

include("gizmos.jl")

k = AbstractAlgebra.GF(2^30+3)

N, D = 6, 6
random = [
    ("rand-$(n)-$d", randsys(n, d))
    for n in 1:N, d in 1:D
]

systems = vcat(random)

println("Running the following systems: ", map(first, systems))

timers = []
for (name, sys) in systems
    @info "Running $name.."
    for t in 1:2
        TimerOutputs.enable_timer!(Groebner._TIMER); 
        reset_timer!(Groebner._TIMER); 
        groebner(sys); 
        if t == 2 show(Groebner._TIMER, allocations=false); println() end
    end
    push!(timers, [name, copy(Groebner._TIMER)])
end

data = Matrix{Any}(undef, N, D)
for i in 1:N, j in 1:D
    name = "rand-$(i)-$j"
    timer = timers[findfirst(t -> t[1] == name, timers)][2]
    time = TimerOutputs.tottime(timer)
    data[i, j] = time
end

pretty_table(
    data,
    column_labels = ["d=$j" for j in 1:D],
    row_labels = ["n=$i" for i in 1:N],
    alignment = :l,
    formatters=[
        (v,i,j) -> v isa Number ? @sprintf("%.2f", v / 1e9) : v, 
    ],
    title = "Time (s) for random systems over $k",
)
