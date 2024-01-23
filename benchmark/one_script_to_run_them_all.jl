# Activate the environment and add all required packages
import Pkg
include("generate/utils.jl")
julia_pkg_preamble("$(@__DIR__)")

# Load the packages
using AbstractAlgebra
using ArgParse
using Base.Threads
using CpuId, Logging, Pkg, Printf
using Distributed
using Dates
using Groebner
using ProgressMeter, PrettyTables
using NaturalSort

# Load the definitions of benchmark systems
include("benchmark_systems.jl")

# Load the code to generate benchmarks for different software
include("generate/benchmark_generators.jl")

# Load the code to compute the certificate of a groebner basis
include("generate/basis_certificate.jl")

# Set up the logger
global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

# Set the properties of progress bars
const _progressbar_color = :light_green
const _progressbar_value_color = :light_green
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

const _available_backends = ["groebner", "singular", "maple", "msolve", "openf4"]

const _skip_singular = true
const _skip_openf4 = true

# Parses command-line arguments
#! format: off
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "backend"
            help = """
            The backend to benchmark.
            Possible options are:
            - groebner
            - learn_apply
            - singular
            - maple
            - msolve
            - openf4
            - ALL"""
            arg_type = String
            required = true
        "--lopenf4"
            help = """
            Required if openf4 is the backend. 
            Must point to the location where openf4 library is installed."""
            arg_type = String
            default = ""
            required = false
        "--bmaple"
            help = """
            If maple is the backend, specifies the location of the maple binary.
            """
            arg_type = String
            default = "maple"
            required = false
        "--bmsolve"
            help = """
            If msolve is the backend, specifies the location of the msolve binary.
            """
            arg_type = String
            default = "msolve"
            required = false
        "--suite"
            help = """
            The index of the benchmark suite.
            Possible options are:
            - 1: for benchmarks over integers modulo a prime
            - 2: for benchmarks over integers modulo a small prime
            - 3: for benchmarks over the rationals"""
            arg_type = Int
            default = 1
        "--validate"
            help = """
            Validate the output bases against the correct ones.
            This will result in a huge slowdown for some of the backends.
                
            Possible options are:
            - `yes`: validate results
            - `no`: do not validate results
            - `update`: update the validation certificates"""
            arg_type = String
            default = "no"
        "--nruns"
            help = "Number of times to run the benchmark."
            arg_type = Int
            default = 1
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 60
        "--nworkers"
            help = "The number of worker processes."
            arg_type = Int
            default = 4
    end

    parse_args(s)
end
#! format: on

function generate_benchmark_file(
    backend,
    name,
    system,
    dir,
    validate,
    nruns,
    time_filename,
    args
)
    if backend == "groebner"
        benchmark_source = generate_benchmark_source_for_groebner(
            name,
            system,
            dir,
            validate,
            nruns,
            time_filename
        )
        fd = open("$dir/$name.jl", "w")
        println(fd, benchmark_source)
        close(fd)
    elseif backend == "learn_apply"
        benchmark_source = generate_benchmark_source_for_groebner(
            name,
            system,
            dir,
            validate,
            nruns,
            time_filename
        )
        fd = open("$dir/$name.jl", "w")
        println(fd, benchmark_source)
        close(fd)
    elseif backend == "singular"
        benchmark_source = generate_benchmark_source_for_singular(
            name,
            system,
            dir,
            validate,
            nruns,
            time_filename
        )
        fd = open("$dir/$name.jl", "w")
        println(fd, benchmark_source)
        close(fd)
    elseif backend == "maple"
        benchmark_source = generate_benchmark_source_for_maple(
            name,
            system,
            dir,
            validate,
            nruns,
            time_filename
        )
        fd = open("$dir/$name.mpl", "w")
        println(fd, benchmark_source)
        close(fd)
    elseif backend == "msolve"
        benchmark_source = generate_benchmark_source_for_msolve(
            name,
            system,
            dir,
            validate,
            nruns,
            time_filename
        )
        fd = open("$dir/$name.in", "w")
        println(fd, benchmark_source)
        close(fd)
    elseif backend == "openf4"
        benchmark_source = generate_benchmark_source_for_openf4(
            name,
            system,
            dir,
            validate,
            nruns,
            time_filename
        )
        fd = open("$dir/$name.cpp", "w")
        println(fd, benchmark_source)
        close(fd)
    end
end

function get_command_to_run_benchmark(
    backend,
    problem_name,
    problem_num_runs,
    problem_set_id,
    validate,
    args
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
    elseif backend == "learn_apply"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/learn_apply/run_in_learn_apply.jl",
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
        return Cmd([args["bmaple"], "$scriptpath/$problem_name/$(problem_name).mpl"])
    elseif backend == "msolve"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/msolve/run_in_msolve.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id",
            "$validate",
            "$(args["bmsolve"])"
        ])
    elseif backend == "openf4"
        return Cmd([
            "julia",
            (@__DIR__) * "/generate/openf4/run_in_openf4.jl",
            "$problem_name",
            "$problem_num_runs",
            "$problem_set_id",
            "$(args["lopenf4"])"
        ])
    end
end

function populate_benchmarks(args; regenerate=true)
    backend = args["backend"]
    benchmark_id = args["suite"]
    nruns = args["nruns"]
    validate = args["validate"] in ["yes", "update"]
    benchmark = get_benchmark_suite(benchmark_id)
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
            time_filename,
            args
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
    benchmark_id = args["suite"]

    benchmark = get_benchmark_suite(benchmark_id)
    benchmark_name = benchmark.name

    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    systems_to_benchmark = first(walkdir(benchmark_dir))[2]
    indices_to_benchmark = collect(1:length(systems_to_benchmark))

    added_timeout_time = 10
    @info """
    Benchmarking $backend.
    Benchmark suite: $benchmark_name
    Number of benchmark systems: $(length(systems_to_benchmark))
    Validate result: $(validate)
    Workers: $(nworkers)
    Timeout: $timeout + $added_timeout_time seconds
    Aggregate results over $nruns runs"""
    @info """
    Benchmark systems:
    $systems_to_benchmark"""

    timeout = timeout + 10  # add 10 seconds for compilation time :)
    seconds_passed(from_t) = round((time_ns() - from_t) / 1e9, digits=2)

    queue = [(problem_id=problem,) for problem in indices_to_benchmark]
    processes = []
    running = []
    errored = []
    timedout = []

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
            cmd = get_command_to_run_benchmark(
                backend,
                problem_name,
                nruns,
                benchmark_id,
                validate,
                args
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
                    kill(proc.julia_process, Base.SIGKILL)
                    close(proc.logfile)
                    # close(proc.errfile)
                    push!(timedout, proc)
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

    if !isempty(timedout)
        printstyled("(!) Timed-out:\n", color=:light_yellow)
        for proc in timedout
            print("$(proc.problem_name), ")
        end
        println()
    end

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

    benchmark_id = args["suite"]
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    validate_dir = (@__DIR__) * "/" * get_validate_dir(benchmark_id)

    @info """Validating results for $backend. May take some time.
    Directory with the certificates is $validate_dir
    """

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
            printstyled("\tMISSING CERTIFICATE...", color=:light_yellow)
        end
        # At this point, the recently computed basis is stored in `result`
        @assert result_exists
        success, result_validation_certificate =
            compute_basis_validation_certificate(result)
        if !success
            @warn "Bad file encountered at $problem_result_path. Skipping"
            continue
        end
        if update_certificates || !true_result_exists
            mkpath("$validate_dir/$problem_name/")
            true_result_file = open(problem_validate_path, "w")
            println(true_result_file, result_validation_certificate)
            printstyled("\tUPDATED\n", color=:light_green)
            continue
        end
        @assert result_exists && true_result_exists
        @assert is_certificate_standardized(result_validation_certificate)
        @assert is_certificate_standardized(true_result)
        if result_validation_certificate != true_result
            printstyled("\tWRONG\n", color=:light_red)
            println("True certificate:\n$true_result")
            println("Current certificate:\n$result_validation_certificate")
        else
            printstyled("\tOK\n", color=:light_green)
        end
    end

    nothing
end

function collect_timings(args, names)
    backend = args["backend"]
    benchmark_id = args["suite"]
    benchmark_dir = (@__DIR__) * "/" * get_benchmark_dir(backend, benchmark_id)
    benchmark_name = get_benchmark_suite(benchmark_id).name

    if backend == "learn_apply"
        targets =
            [:total_time_F4, :total_time_learn, :total_time_apply, :total_time_apply_4x]
    else
        targets = [:total_time]
    end
    @assert length(targets) > 0
    println()
    @info """
    Collecting results for $backend. Statistics of interest:
    \t$(join(map(string, targets), "\n\t"))
    """

    cannot_collect = []
    names = sort(names, lt=NaturalSort.natural)

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
            continue
        end
        lines = readlines(data_file)
        if isempty(lines)
            @debug "Cannot collect data for $name"
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
        printstyled("(!) Cannot collect benchmark data for:\n", color=:light_yellow)
        for (name,) in cannot_collect
            print("$name, ")
        end
        println()
    end

    title = "Benchmark results, $backend"
    conf = set_pt_conf(tf=tf_markdown, alignment=:c)
    makecolname(target) = HUMAN_READABLE_CATEGORIES[target]
    columns = [makecolname(target) for target in targets]
    header = vcat("Model", map(string, columns))
    vec_of_vecs = Vector{Vector{Any}}()
    for name in names
        push!(vec_of_vecs, [name])
        for target in targets
            formatting_style = CATEGORY_FORMAT[target]
            if haskey(runtime, name) && haskey(runtime[name], target)
                push!(vec_of_vecs[end], formatting_style(runtime[name][target]))
            else
                push!(vec_of_vecs[end], "-")
            end
        end
    end

    table = Array{Any, 2}(undef, length(vec_of_vecs), length(header))
    for i in 1:length(vec_of_vecs)
        for j in 1:length(vec_of_vecs[i])
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

    columns = [makecolname(target) for target in targets]
    resulting_md *= "|Model|" * join(map(string, columns), "|") * "|\n"
    resulting_md *= "|:----|" * join(["---" for _ in columns], "|") * "|\n"
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

    return runtime
end

function collect_all_timings(args, runtimes, systems)
    benchmark_id = args["suite"]
    benchmark_name = get_benchmark_suite(benchmark_id).name

    targets = [:total_time]
    @assert length(targets) > 0

    backends = collect(keys(runtimes))
    backends = sort(backends, lt=NaturalSort.natural)
    systems = sort(systems, lt=NaturalSort.natural)

    println()
    @info """
    Collecting all results. Statistics of interest:
    \t$(join(map(string, targets), "\n\t"))
    Present backends:
    \t$(join(map(string, backends), "\n\t"))
    """

    _target = targets[1]
    formatting_style = CATEGORY_FORMAT[_target]
    conf = set_pt_conf(tf=tf_markdown, alignment=:c)
    title = "Benchmark results"
    header = vcat("System", backends)
    vec_of_vecs = Vector{Vector{Any}}()
    for name in systems
        row = [name]
        for backend in backends
            if haskey(runtimes[backend], name) && haskey(runtimes[backend][name], _target)
                push!(row, formatting_style(runtimes[backend][name][_target]))
            else
                push!(row, "-")
            end
        end
        push!(vec_of_vecs, row)
    end
    table = Array{Any, 2}(undef, length(vec_of_vecs), length(backends) + 1)
    for i in 1:length(vec_of_vecs)
        for j in 1:(length(backends) + 1)
            table[i, j] = vec_of_vecs[i][j]
        end
    end
    println()
    h1 = Highlighter(
        (data, i, j) ->
            j > 1 &&
                strip(data[i, j]) != "-" &&
                parse(Float64, data[i, j]) == minimum(
                    map(
                        x -> parse(Float64, x),
                        filter(x -> strip(x) != "-", data[i, 2:end])
                    )
                ),
        bold=true
    )
    pretty_table_with_conf(
        conf,
        table;
        header=header,
        title=title,
        limit_printing=false,
        highlighters=(h1,)
    )
    println("All results are in seconds.")

    # Print the table to BENCHMARK_TABLE.
    resulting_md = ""
    resulting_md *= """
    ## Benchmark results

    $(now())

    Benchmarked backends: $backends

    Benchmark suite: $benchmark_name

    - Workers: $(args["nworkers"])
    - Timeout: $(args["timeout"]) s
    - Aggregated over: $(args["nruns"]) runs

    **All timings in seconds.**

    """

    makecolname(target) = HUMAN_READABLE_CATEGORIES[target]
    columns = backends
    resulting_md *= "|Model|" * join(map(string, columns), "|") * "|\n"
    resulting_md *= "|:----|" * join(["---" for _ in columns], "|") * "|\n"
    for name in systems
        resulting_md *= "|$name|"
        for backend in backends
            model_data = runtimes[backend][name]
            for target in targets
                if !haskey(model_data, target)
                    resulting_md *= " - " * "|"
                else
                    formatting_style = CATEGORY_FORMAT[target]
                    resulting_md *= formatting_style(model_data[target]) * "|"
                end
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
        (@__DIR__) * "/$BENCHMARK_RESULTS/$(BENCHMARK_TABLE)_$(benchmark_id).md"
    open(table_filename, "w") do io
        write(io, resulting_md)
    end

    println()
    printstyled("Table with results", color=:light_green)
    println(" is written to $table_filename")
end

function check_args(args)
    backend = args["backend"]
    @assert backend in ("groebner", "singular", "maple", "openf4", "msolve", "ALL")
    if backend == "openf4" && args["suite"] in [3, 7]
        throw("Running benchmarks over the rationals is not possible for openf4")
    end
    # if backend == "msolve" && args["suite"] in [3, 7] && args["validate"] != "no"
    #     throw(
    #         "Validating results for msolve over the rationals is not possible. Use command line option --validate=no"
    #     )
    # end
    if backend == "learn_apply" && args in [3, 7]
        throw("Cannot learn & apply over the rationals")
    end
end

function cleanup(names=nothing)
    # this is quite bad but fine :^)
    cmd1 = Cmd(`pkill msolve`, ignorestatus=true)
    cmd2 = Cmd(`pkill mserver`, ignorestatus=true)
    run(cmd1)
    run(cmd2)
    if !(names === nothing)
        for name in names
            cmd = Cmd(`pkill -f "$name"`, ignorestatus=true)
            run(cmd)
        end
    end
    nothing
end

function main()
    # Parse command line args
    args = parse_commandline()
    @debug "Command-line args:"
    for (arg, val) in args
        @debug "$arg  =>  $val"
    end

    println("Timestamp: $(now())")

    # Either benchmark all available backends or a single backend
    if args["backend"] == "ALL"
        @info "Benchmarking all available backends"
        runtimes = Dict()
        systems = []
        solved_problems = []
        for backend in _available_backends
            if _skip_singular && backend == "singular"
                continue
            end
	    if _skip_openf4 && backend == "openf4"
                continue
            end
            args_ = copy(args)
            args_["backend"] = backend
            try
                check_args(args_)
                populate_benchmarks(args_)
                solved_problems = run_benchmarks(args_)
                validate_results(args_, solved_problems)
                runtime = collect_timings(args_, solved_problems)
                runtimes[backend] = runtime
                union!(systems, solved_problems)
            catch e
                printstyled("(!) ", color=:light_yellow)
                println("Cannot benchmark $backend")
            finally
                cleanup(systems)
            end
        end
        collect_all_timings(args, runtimes, systems)
    else
        solved_problems = []
        try
            check_args(args)
            populate_benchmarks(args)
            solved_problems = run_benchmarks(args)
            validate_results(args, solved_problems)
            collect_timings(args, solved_problems)
        finally
            cleanup(solved_problems)
        end
    end
end

main()
