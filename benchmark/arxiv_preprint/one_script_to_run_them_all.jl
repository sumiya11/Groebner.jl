# Load all required packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using ArgParse
using CpuId, Logging, Pkg, Printf
using Base.Threads
using Distributed
using Dates
using ProgressMeter

using Groebner, AbstractAlgebra

global_logger(Logging.ConsoleLogger(stdout, Logging.Debug))

include("benchmark_systems.jl")
include("utils.jl")

const _progressbar_color = :light_green
const _progressbar_value_color = :light_green
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

const BENCHMARK_TABLE = "benchmark_result"

function parse_commandline()
    s = ArgParseSettings()
    #! format: off
    @add_arg_table s begin
        "backend"
            help = """
            The Groebner basis computation backend to benchmark.
            Possible options are:
            - groebner
            - singular
            - maple
            - msolve
            - openf4"""
            arg_type = String
            required = true
        "--benchmark"
            help = """
            The index of benchmark dataset.
            Possible options are:
            - 1: for benchmarks over integers modulo 2^31-1
            - 2: for benchmarks over the rationals"""
            arg_type = Int
            default = 1
        "--nruns"
            help = "Number of times to run the benchmark."
            arg_type = Int
            default = 1
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 300
        "--nworkers"
            help = "The number of worker processes."
            arg_type = Int
            default = 4
    end
    #! format: on

    parse_args(s)
end

function generate_benchmark_file(backend, name, system, dir)
    if backend == "groebner"
        # generate strings to write to a file
        ring = parent(system[1])
        field = base_ring(ring)
        vars_repr = join(map(string, gens(ring)), ", ")
        field_repr = if iszero(characteristic(field))
            "QQ"
        else
            "GF($(characteristic(field)))"
        end
        vars_repr_quoted = map(s -> "$s", map(string, gens(ring)))
        ring_repr = """ring, ($vars_repr) = PolynomialRing(
            $field_repr, 
            $vars_repr_quoted, 
            ordering=:degrevlex
        )"""
        system_repr = join(map(s -> "\t" * s, map(repr, system)), ",\n")
        system_repr = "system = [\n$system_repr\n]"
        fd = open("$dir/$name.jl", "w")
        println(fd, "# $name")
        println(fd, "#! format: off")
        println(fd, "using AbstractAlgebra, Groebner")
        println(fd, "")
        println(fd, ring_repr)
        println(fd, system_repr)
        close(fd)
    end
end

function command_to_run_a_single_system(
    backend,
    problem_name,
    problem_num_runs,
    problem_set_id
)
    if backend == "groebner"
        return Cmd([
            "julia",
            (@__DIR__) * "/run_in_julia.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id"
        ])
    end
end

function populate_benchmarks(args; regenerate=true)
    backend = args["backend"]
    benchmark_id = args["benchmark"]
    benchmark = get_benchmark(benchmark_id)
    benchmark_name, systems = benchmark.name, benchmark.systems
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(benchmark_id)
    dir_present = isdir(benchmark_dir)
    if !dir_present || regenerate
        @info "Re-generating the folder with benchmarks"
        try
            if isdir(benchmark_dir)
                rm(benchmark_dir, recursive=true, force=true)
            end
        catch err
            @info "Something went wrong when deleting the benchmarks folder"
            showerror(stdout, err)
            println(stdout)
        end
    end
    prog = Progress(
        length(systems),
        "Generating benchmark files: $benchmark_name",
        # spinner = true,
        dt=0.1,
        enabled=progressbar_enabled(),
        color=_progressbar_color
    )
    for bmark in systems
        next!(prog) # , spinner = "⌜⌝⌟⌞")
        system_name = bmark[1]
        system = bmark[2]
        @debug "Generating $system_name"
        benchmark_system_dir = "$benchmark_dir/$system_name/"
        mkpath(benchmark_system_dir)
        generate_benchmark_file(backend, system_name, system, benchmark_system_dir)
    end
    finish!(prog)
    true
end

function run_benchmarks(args)
    timestamp = time_ns()
    timeout = args["timeout"]
    @assert timeout > 0
    backend = args["backend"]
    @assert backend in ("groebner", "singular", "maple", "openf4", "msolve")
    nworkers = args["nworkers"]
    @assert nworkers > 0
    nruns = args["nruns"]
    @assert nruns > 0
    benchmark_id = args["benchmark"]

    benchmark = get_benchmark(benchmark_id)
    benchmark_name = benchmark.name

    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(benchmark_id)
    systems_to_benchmark = first(walkdir(benchmark_dir))[2]
    indices_to_benchmark = collect(1:length(systems_to_benchmark))

    @info """
    Benchmarking $backend."""
    @info """
    Number of benchmark systems: $(length(systems_to_benchmark))
    Workers: $(nworkers)
    Timeout: $timeout seconds"""
    @info """
    Benchmark systems:
    $systems_to_benchmark"""

    seconds_passed(from_t) = round((time_ns() - from_t) / 1e9, digits=2)

    queue = [(problem_id=problem,) for problem in indices_to_benchmark]
    processes = []
    running = []
    errored = []

    generate_showvalues(processes) =
        () -> [(
            :Active,
            join(
                map(
                    proc -> string(proc.problem_name),
                    filter(proc -> process_running(proc.julia_process), processes)
                ),
                ", "
            )
        )]

    prog = Progress(
        length(queue),
        "Running benchmarks",
        # spinner = true,
        dt=0.3,
        enabled=progressbar_enabled(),
        color=_progressbar_color
    )
    while true
        if !isempty(queue) && length(running) < nworkers
            task = pop!(queue)
            problem_id = task.problem_id
            problem_name = systems_to_benchmark[problem_id]
            log_filename = generic_filename("logs")
            log_file = open("$benchmark_dir/$problem_name/$log_filename", "w")
            @debug "Running $problem_name. Logs: $benchmark_dir/$problem_name/$log_filename"
            cmd = command_to_run_a_single_system(backend, problem_name, nruns, benchmark_id)
            cmd = Cmd(cmd, ignorestatus=true, detach=false, env=copy(ENV))
            proc = run(pipeline(cmd, stdout=log_file, stderr=log_file), wait=false)
            push!(
                processes,
                (
                    problem_id=problem_id,
                    problem_name=problem_name,
                    julia_process=proc,
                    start_time=time_ns(),
                    logfile=log_file
                )
            )
            push!(running, processes[end])
            next!(
                prog,
                showvalues = generate_showvalues(processes),
                step       = 0,
                valuecolor = _progressbar_value_color
                # spinner = "⌜⌝⌟⌞",
            )
        end

        sleep(0.2)
        to_be_removed = []
        for i in 1:length(running)
            proc = running[i]
            if process_exited(proc.julia_process)
                push!(to_be_removed, i)
                if proc.julia_process.exitcode != 0
                    push!(errored, proc)
                end
                close(proc.logfile)
                # close(proc.errfile)
                start_time = proc.start_time
                next!(
                    prog,
                    showvalues = generate_showvalues(processes),
                    valuecolor = _progressbar_value_color
                )
                @debug "Yielded $(proc.problem_name) after $(seconds_passed(start_time)) seconds"
            end
            if process_running(proc.julia_process)
                start_time = proc.start_time
                if seconds_passed(start_time) > timeout
                    push!(to_be_removed, i)
                    kill(proc.julia_process)
                    close(proc.logfile)
                    # close(proc.errfile)
                    next!(
                        prog,
                        showvalues = generate_showvalues(processes),
                        valuecolor = _progressbar_value_color
                    )
                    @debug "Timed-out $(proc.problem_name) after $(seconds_passed(start_time)) seconds"
                end
            end
        end
        deleteat!(running, to_be_removed)
        if isempty(queue) && isempty(running)
            @debug "All benchmarks finished"
            break
        end
    end
    finish!(prog)

    if !isempty(errored)
        printstyled("(!) Maybe errored:\n", color=:red)
        for proc in errored
            println("\t$(proc.problem_name)")
        end
    end

    printstyled(
        """
        Benchmarking finished in $(round((time_ns() - timestamp) / 1e9, digits=2)) seconds.
        Results are written to $benchmark_dir
        """,
        color=:light_green
    )

    systems_to_benchmark
end

function collect_timings(args, names; content=:compare)
    backend = args["backend"]
    benchmark_id = args["benchmark"]
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(benchmark_id)

    targets = [:total_time]
    @assert length(targets) > 0
    @info """
    Collecting benchmark results for $backend.

    Statistics of interest:
    \t$(join(map(string, targets), "\n\t"))
    """

    cannot_collect = []
    names = sort(names)

    # Collect timings and data from directory BENCHMARK_RESULTS.
    runtime = Dict()
    for name in names
        @debug "==== Reading $name"
        runtime[name] = Dict()
        timingsfn = timings_filename()
        timings_file = nothing
        #####
        try
            @debug "==== Opening $benchmark_dir/$name/$timingsfn"
            timings_file = open("$benchmark_dir/$name/$timingsfn", "r")
        catch e
            @debug "Cannot collect timings for $name"
            push!(cannot_collect, (name,))
            continue
        end
        lines = readlines(timings_file)
        if isempty(lines)
            @debug "Cannot collect timings for $name"
            push!(cannot_collect, (name,))
            continue
        end
        @assert lines[1] == name
        for line in lines[2:end]
            k, v = split(line, ", ")
            runtime[name][Symbol(k)] = parse(Float64, v)
        end
        close(timings_file)
        #####
        datafn = generic_filename("data")
        data_file = nothing
        try
            @debug "==== Opening $benchmark_dir/$name/$datafn"
            data_file = open("$benchmark_dir/$name/$datafn", "r")
        catch e
            @debug "Cannot collect data for $name"
            push!(cannot_collect, (name,))
            continue
        end
        lines = readlines(data_file)
        if isempty(lines)
            @debug "Cannot collect data for $name"
            push!(cannot_collect, (name,))
            continue
        end
        @assert lines[1] == name
        for line in lines[2:end]
            k, v = map(strip, split(line, ","))
            runtime[name][Symbol(k)] = v
        end
        close(data_file)
    end

    if !isempty(cannot_collect)
        printstyled("(!) Cannot collect benchmark data for:\n", color=:red)
        for (name,) in cannot_collect
            println("\t$name")
        end
    end

    # Print the table to BENCHMARK_TABLE.
    resulting_md = ""
    resulting_md *= """
    ## Benchmark results

    $(now())

    - Benchmarked backend: `$backend`
    - Workers: $(args["nworkers"])
    - Timeout: $(args["timeout"]) s

    **All timings in seconds.**

    """

    makecolname(target) = HUMAN_READABLE_CATEGORIES[target]
    columns = [makecolname(target) for target in targets]
    resulting_md *= "|Model|" * join(map(string, columns), "|") * "|\n"
    resulting_md *= "|-----|" * join(["---" for _ in columns], "|") * "|\n"
    for name in names
        model_data = runtime[name]
        resulting_md *= "|$name|"
        for target in targets
            if !haskey(model_data, target)
                resulting_md *= " - " * "|"
            else
                formatting_style = CATEGORY_FORMAT[target]
                resulting_md *= formatting_style(model_data[target]) * "|"
            end
        end
        resulting_md *= "\n"
    end

    resulting_md *= "\n*Benchmarking environment:*\n\n"
    resulting_md *= "* Total RAM (GiB): $(div(Sys.total_memory(), 2^30))\n"
    resulting_md *= "* Processor: $(cpubrand())\n"
    resulting_md *= "* Julia version: $(VERSION)\n\n"
    resulting_md *= "Versions of the dependencies:\n\n"

    deps = Pkg.dependencies()
    stid_info = deps[findfirst(x -> x.name == "StructuralIdentifiability", deps)]
    for (s, uid) in stid_info.dependencies
        if deps[uid].version !== nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    table_filename =
        (@__DIR__) * "/$BENCHMARK_RESULTS/$(BENCHMARK_TABLE)_$(benchmark_id).md"
    open(table_filename, "w") do io
        write(io, resulting_md)
    end

    printstyled("Table with results is written to $table_filename", color=:light_green)
end

function main()
    # Parse command line args
    args = parse_commandline()
    @debug "Command-line args:"
    for (arg, val) in args
        @debug "$arg  =>  $val"
    end

    # Create directories with benchmarks
    populate_benchmarks(args)

    # Run benchmarks and record results
    problems = run_benchmarks(args)

    # Collect the timings and other info
    collect_timings(args, problems)
end

main()