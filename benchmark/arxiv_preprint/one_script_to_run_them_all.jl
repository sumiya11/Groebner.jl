# Add all required packages
import Pkg

include("generate/utils.jl")
julia_pkg_preamble(".")

# Load the packages
using ArgParse
using Base.Threads
using CpuId, Logging, Pkg, Printf
using Distributed
using Dates
using ProgressMeter, PrettyTables
using AbstractAlgebra, Groebner

# Set the logger
global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

# Load benchmark systems
include("generate/benchmark_systems.jl")

# Load the code to compute the certificate of a groebner basis
include("generate/basis_certificate.jl")

# Set the properties of progress bar
const _progressbar_color = :light_green
const _progressbar_value_color = :light_green
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

# Parses command-line arguments
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
        "lopenf4"
            help = """
            Required if openf4 is specified as the backend. 
            Must point to the location where openf4 library is installed."""
            arg_type = String
            default = ""
            required = false
        "--benchmark"
            help = """
            The index of benchmark dataset.
            Possible options are:
            - 1: for benchmarks over integers modulo a prime
            - 2: for benchmarks over the rationals"""
            arg_type = Int
            default = 1
        "--validate"
            help = """
            Validate the output bases against the (presumably) correct ones.
            This can result in a slowdown for some of the backends.
                
            Possible options are:
            - `yes`: validate results
            - `no`: do not validate results
            - `update`: update validation certificates"""
            arg_type = String
            default = "yes"
        "--nruns"
            help = "Number of times to run the benchmark."
            arg_type = Int
            default = 1
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 600
        "--nworkers"
            help = "The number of worker processes."
            arg_type = Int
            default = 4
    end
    #! format: on

    parse_args(s)
end

function generate_benchmark_file(backend, name, system, dir, validate, nruns, time_filename)
    ring = parent(system[1])
    field = base_ring(ring)
    if backend == "groebner" || backend == "singular"
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
    elseif backend == "maple"
        fd = open("$dir/$name.mpl", "w")
        println(fd, "# $name")
        println(fd, "with(Groebner):")
        println(fd, "with(PolynomialIdeals):")
        println(fd, "kernelopts(numcpus=1);")
        system_repr = join(map(s -> "\t\t" * s, map(repr, system)), ",\n")
        vars_repr = join(map(string, gens(ring)), ", ")
        println(fd, "")
        println(fd, "runtime := 2^1000:")
        println(fd, "for i from 1 by 1 to $nruns do")
        println(fd, "\tJ := [\n$system_repr\n\t]:")
        println(fd, "\tprint(\"Running $name\");")
        println(fd, "\tst := time[real]():")
        println(
            fd,
            "\tG := Groebner[Basis](J, tdeg($vars_repr), method=fgb, characteristic=$(characteristic(field))):"
        )
        println(fd, "\tprint(\"$name: \", time[real]() - st);")
        println(fd, "\truntime := min(runtime, time[real]() - st);")
        println(fd, "end do;")
        println(fd, "")
        println(fd, "timings_fn := \"$time_filename\":")
        println(fd, "FileTools[Text][WriteLine](timings_fn, \"$name\");")
        println(
            fd,
            "FileTools[Text][WriteLine](timings_fn, cat(\"total_time, \", String(runtime)));"
        )
        if validate
            println(fd)
            output_fn = output_filename()
            println(fd, "output_fn := \"$dir/$output_fn\":")
            println(fd, "FileTools[Text][WriteLine](output_fn, \"$vars_repr\");")
            println(
                fd,
                "FileTools[Text][WriteLine](output_fn, \"$(characteristic(field))\");"
            )
            println(
                fd,
                """
                for poly in G do
                    FileTools[Text][WriteLine](output_fn, cat(String(poly), \",\")):
                end do;
                """
            )
        end
        close(fd)
    elseif backend == "msolve"
        vars_repr = join(map(string, gens(ring)), ", ")
        system_repr = join(map(repr, system), ",\n")
        fd = open("$dir/$name.in", "w")
        println(fd, "$vars_repr")
        println(fd, "$(characteristic(field))")
        println(fd, system_repr)
        close(fd)
    elseif backend == "openf4"
        fd = open("$dir/$name.cpp", "w")
        println(fd, "// $name")
        println(fd, "#include <iostream>")
        println(fd, "#include <fstream>")
        println(fd, "#include <string>")
        println(fd, "#include <vector>")
        println(fd, "#include <libopenf4.h>")
        println(
            fd,
            """
            using namespace std;

            int main (int argc, char **argv)
            {
            """
        )
        vars_repr =
            join(map(s -> "\tvariableName.push_back(\"$s\");", map(repr, gens(ring))), "\n")
        system_repr = join(
            map(s -> "\tpolynomialArray.emplace_back(\"$s\");", map(repr, system)),
            "\n"
        )
        println(
            fd,
            """
            \tvector<string> polynomialArray;
            \tvector<string> variableName;

            $vars_repr
            $system_repr

            \tvector<string> basis = groebnerBasisF4($(characteristic(field)), $(length(gens(ring))), variableName, polynomialArray, 1, 0);

            \tstd::cout << \"The basis contains \" << basis.size() << \" elements.\" << std::endl;
            """
        )
        if validate
            vars_ = join(map(repr, gens(ring)), ", ")
            output_fn = output_filename()
            println(fd)
            println(
                fd,
                """
                \tofstream output;
                \toutput.open(\"$dir/$output_fn\");
                \toutput << \"$vars_\" << endl;
                \toutput << \"$(characteristic(field))\" << endl;
                \tfor (size_t i = 0; i < basis.size(); i++) {
                \t\toutput << basis[i] << \",\" << endl;
                \t}
                \toutput.close();
                """
            )
        end
        println(fd, "\treturn 0;")
        println(fd, "}")
        close(fd)
    end
end

function command_to_run_a_single_system(
    backend,
    problem_name,
    problem_num_runs,
    problem_set_id,
    validate;
    lib=nothing
)
    if backend == "groebner"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/groebner/run_in_groebner.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id",
            "$validate"
        ])
    elseif backend == "singular"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/singular/run_in_singular.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id",
            "$validate"
        ])
    elseif backend == "maple"
        scriptpath = (@__DIR__) * "/" * get_benchmark_dir(backend, problem_set_id)
        return Cmd(["maple", "$scriptpath/$problem_name/$(problem_name).mpl"])
    elseif backend == "msolve"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/msolve/run_in_msolve.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id",
            "$validate"
        ])
    elseif backend == "openf4"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/openf4/run_in_openf4.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id",
            "$lib"
        ])
    end
    throw("Unknown backend")
end

function populate_benchmarks(args; regenerate=true)
    backend = args["backend"]
    benchmark_id = args["benchmark"]
    nruns = args["nruns"]
    validate = args["validate"] in ["yes", "update"]
    benchmark = get_benchmark(benchmark_id)
    benchmark_name, systems = benchmark.name, benchmark.systems
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    dir_present = isdir(benchmark_dir)
    if !dir_present || regenerate
        @info "Re-generating the folder with benchmarks"
        try
            if isdir(benchmark_dir)
                rm(benchmark_dir, recursive=true, force=true)
            end
        catch err
            @info "Something went wrong when deleting the directory with benchmarks"
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
        time_filename = "$benchmark_dir/$system_name/$(timings_filename())"
        generate_benchmark_file(
            backend,
            system_name,
            system,
            benchmark_system_dir,
            validate,
            nruns,
            time_filename
        )
    end
    finish!(prog)
    true
end

function run_benchmarks(args)
    timestamp = time_ns()
    timeout = args["timeout"]
    @assert timeout > 0
    backend = args["backend"]
    nworkers = args["nworkers"]
    validate = args["validate"] in ["yes", "update"]
    @assert nworkers > 0
    nruns = args["nruns"]
    @assert nruns > 0
    benchmark_id = args["benchmark"]

    benchmark = get_benchmark(benchmark_id)
    benchmark_name = benchmark.name

    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    systems_to_benchmark = first(walkdir(benchmark_dir))[2]
    indices_to_benchmark = collect(1:length(systems_to_benchmark))

    @info """
    Benchmarking $backend.
    Benchmark suite: $benchmark_name
    Number of benchmark systems: $(length(systems_to_benchmark))
    Validate result: $(validate)
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
            cmd = command_to_run_a_single_system(
                backend,
                problem_name,
                nruns,
                benchmark_id,
                validate,
                lib=args["lopenf4"]
            )
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
                showvalues = generate_showvalues(running),
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
                    showvalues = generate_showvalues(running),
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
                        showvalues = generate_showvalues(running),
                        valuecolor = _progressbar_value_color
                    )
                    @debug "Timed-out $(proc.problem_name) after $(seconds_passed(start_time)) seconds"
                end
            end
        end
        deleteat!(running, to_be_removed)
        next!(
            prog,
            showvalues = generate_showvalues(running),
            step       = 0,
            valuecolor = _progressbar_value_color
            # spinner = "⌜⌝⌟⌞",
        )
        if isempty(queue) && isempty(running)
            @debug "All benchmarks finished"
            break
        end
    end
    finish!(prog)

    if !isempty(errored)
        printstyled("(!) Maybe errored:\n", color=:light_red)
        for proc in errored
            print("$(proc.problem_name), ")
        end
        println()
    end

    println()
    println(
        "Benchmarking finished in $(round((time_ns() - timestamp) / 1e9, digits=2)) seconds."
    )
    printstyled("Benchmark results", color=:light_green)
    println(" are written to $benchmark_dir")

    systems_to_benchmark
end

function validate_results(args, problem_names)
    println()
    backend = args["backend"]
    if !(args["validate"] in ["yes", "update"])
        @info "Skipping result validation for $backend"
        return nothing
    end

    update_certificates = args["validate"] == "update"

    benchmark_id = args["benchmark"]
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    validate_dir = (@__DIR__) * "/" * get_validate_dir(benchmark_id)

    @info "Validating results for $backend. May take some time."

    if update_certificates
        @info "Re-generating the folder with certificates"
        try
            if isdir(validate_dir)
                rm(validate_dir, recursive=true, force=true)
            end
        catch err
            @info "Something went wrong when deleting the directory with certificates"
            showerror(stdout, err)
            println(stdout)
        end
    end

    for problem_name in problem_names
        print("$problem_name:")
        problem_validate_path = "$validate_dir/$problem_name/$(certificate_filename())"
        problem_result_path = "$benchmark_dir/$problem_name/$(output_filename())"
        true_result_exists, true_result = false, nothing
        result_exists, result = false, nothing
        try
            result_file = open(problem_result_path, "r")
            result = read(result_file, String)
            result_exists = true
            if isempty(string(strip(result, [' ', '\n', '\r'])))
                result_exists = false
            end
        catch e
            @debug "Cannot collect result data for $name"
        end
        if !result_exists
            printstyled("\tMISSING RESULT\n", color=:light_yellow)
            continue
        end
        try
            true_result_file = open(problem_validate_path, "r")
            true_result = read(true_result_file, String)
            true_result = standardize_certificate(true_result)
            true_result_exists = true
        catch e
            @debug "Cannot collect validation data for $name"
            printstyled("\tMISSING CERTIFICATE.. ", color=:light_yellow)
        end
        # At this point, the recently computed basis is stored in `result`
        @assert result_exists
        result_validation_hash = compute_basis_validation_hash(result)
        if update_certificates || !true_result_exists
            mkpath("$validate_dir/$problem_name/")
            true_result_file = open(problem_validate_path, "w")
            println(true_result_file, result_validation_hash)
            printstyled("\tUPDATED\n", color=:light_green)
            continue
        end
        @assert result_exists && true_result_exists
        @assert is_certificate_standardized(result_validation_hash)
        @assert is_certificate_standardized(true_result)
        if result_validation_hash != true_result
            printstyled("\tWRONG RESULT\n", color=:light_red)
            println("True certificate:\n$true_result")
            println("Current certificate:\n$result_validation_hash")
        else
            printstyled("\tOK\n", color=:light_green)
        end
    end

    nothing
end

function collect_timings(args, names)
    backend = args["backend"]
    benchmark_id = args["benchmark"]
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    benchmark_name = get_benchmark(benchmark_id).name

    targets = [:total_time]
    @assert length(targets) > 0
    println()
    @info """
    Collecting results for $backend.
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
            # push!(cannot_collect, (name,))
            continue
        end
        lines = readlines(data_file)
        if isempty(lines)
            @debug "Cannot collect data for $name"
            # push!(cannot_collect, (name,))
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
        printstyled("(!) Cannot collect benchmark data for:\n", color=:light_red)
        for (name,) in cannot_collect
            print("$name, ")
        end
        println()
    end

    _target = targets[1]
    formatting_style = CATEGORY_FORMAT[_target]
    conf = set_pt_conf(tf=tf_markdown, alignment=:c)
    title = "Benchmark results, $backend"
    header = ["System", "Time, s"]
    vec_of_vecs = Vector{Vector{Any}}()
    for name in names
        if haskey(runtime, name) && haskey(runtime[name], _target)
            push!(vec_of_vecs, [name, formatting_style(runtime[name][_target])])
        else
            push!(vec_of_vecs, [name, "-"])
        end
    end
    table = Array{Any, 2}(undef, length(vec_of_vecs), 2)
    for i in 1:length(vec_of_vecs)
        for j in 1:2
            table[i, j] = vec_of_vecs[i][j]
        end
    end
    println()
    pretty_table_with_conf(conf, table; header=header, title=title, limit_printing=false)

    # Print the table to BENCHMARK_TABLE.
    resulting_md = ""
    resulting_md *= """
    ## Benchmark results

    $(now())

    Benchmarked backend: $backend
    Benchmark suite: $benchmark_name

    - Workers: $(args["nworkers"])
    - Timeout: $(args["timeout"]) s
    - Aggregated over: $(args["nruns"]) runs

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
    stid_info = deps[findfirst(x -> x.name == "Groebner", deps)]
    for (s, uid) in stid_info.dependencies
        if deps[uid].version !== nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    table_filename =
        (@__DIR__) * "/$BENCHMARK_RESULTS/$backend/$(BENCHMARK_TABLE)_$(benchmark_id).md"
    open(table_filename, "w") do io
        write(io, resulting_md)
    end

    println()
    printstyled("Table with results", color=:light_green)
    println(" is written to $table_filename")
end

function check_args(args)
    backend = args["backend"]
    @assert backend in ("groebner", "singular", "maple", "openf4", "msolve")
end

function main()
    # Parse command line args
    args = parse_commandline()
    @debug "Command-line args:"
    for (arg, val) in args
        @debug "$arg  =>  $val"
    end

    check_args(args)

    # Create directories with benchmarks
    populate_benchmarks(args)

    # Run benchmarks and record results
    computed_problems = run_benchmarks(args)

    # Validate computed Groebner bases
    validate_results(args, computed_problems)

    # Collect the timings and other info
    collect_timings(args, computed_problems)
end

main()
