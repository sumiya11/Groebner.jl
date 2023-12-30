# This file is a part of Groebner.jl. License is GNU GPL v2.

# Logging for Groebner

# Provides the macro @log, which can be used to log a record with the given
# message to console. The macro @log wraps Logging.@logmsg from the standard
# library.

# Logging is disabled on all threads except the one with threadid() == 1.

"""
    @log expr
    @log level=N expr

Logs a record with `expr` as a message.
Allows to specify the logging level with `level=N`.

## Examples

```jldoctest
@log "Hello, world!"
@log level=1 "Hello, world!"
```
"""
macro log end

###
# Some aux. definitions

function meta_formatter_groebner end

@static if VERSION >= v"1.7.0"
    backend_logger(io, level) = Logging.ConsoleLogger(
        io,
        level,
        show_limited=false,
        meta_formatter=meta_formatter_groebner
    )
else
    backend_logger(io, level) = Logging.ConsoleLogger(io, level)
end

###
# Groebner logger

struct GroebnerLogger <: Logging.AbstractLogger
    logger::ConsoleLogger

    GroebnerLogger() = GroebnerLogger(stderr, Logging.LogLevel(0))
    function GroebnerLogger(io, loglevel::Logging.LogLevel)
        new(backend_logger(io, loglevel))
    end
end

Logging.min_enabled_level(logger::GroebnerLogger) = Logging.min_enabled_level(logger.logger)

function Logging.shouldlog(logger::GroebnerLogger, level, _module, group, id)
    Logging.shouldlog(logger.logger, level, _module, group, id)
end

Logging.catch_exceptions(logger::GroebnerLogger) = Logging.catch_exceptions(logger.logger)

function Logging.handle_message(
    logger::GroebnerLogger,
    lvl,
    msg,
    _mod,
    group,
    id,
    file,
    line;
    kwargs...
)
    # TODO: consider using ActiveFilteredLogger from LoggingExtras.jl
    Logging.handle_message(logger.logger, lvl, msg, _mod, group, id, file, line; kwargs...)
    nothing
end

function meta_formatter_groebner(level::LogLevel, _module, group, id, file, line)
    @nospecialize
    color = Logging.default_logcolor(level)
    prefix = if level >= Logging.Warn
        "Warning"
    elseif level < Logging.Warn && level >= Logging.Info
        "Info"
    else
        "Debug"
    end
    prefix = string(prefix, ":")
    suffix::String = ""
    Logging.Info <= level < Logging.Warn && return color, prefix, suffix
    _module !== nothing && (suffix *= string(_module)::String)
    if file !== nothing
        _module !== nothing && (suffix *= " ")
        suffix *= Base.contractuser(file)::String
        if line !== nothing
            suffix *= ":$(isa(line, UnitRange) ? "$(first(line))-$(last(line))" : line)"
        end
    end
    !isempty(suffix) && (suffix = "@ " * suffix)
    return color, prefix, suffix
end

const _groebner_logger = Ref{GroebnerLogger}(GroebnerLogger())

# Updates the global logging parameters in the Groebner module. 
function update_logger(; loglevel=nothing)
    # Do nothing if logging is disabled
    !logging_enabled() && return nothing
    # Do nothing if run from a worker thread
    threadid() != 1 && return nothing
    if loglevel !== nothing
        _groebner_logger[] = Groebner.GroebnerLogger(stderr, Logging.LogLevel(loglevel))
    end
    nothing
end

###
# The @log macro

@noinline __throw_log_macro_error(file, line, error) =
    throw(ArgumentError("Invalid syntax for @log macro. $error"))

const _default_message_loglevel = Logging.Info

# Parses the arguments to the @log macro and returns a tuple of expressions that
# would evaluate to (loglevel, msg1, msg2, ...).
function log_macro_pruneargs(file, line, args)
    length(args) < 1 &&
        __throw_log_macro_error(file, line, "Argument list must be non-empty")
    level = _default_message_loglevel
    length(args) == 1 && return level, args
    if args[1] isa Expr && args[1].head == :(=)
        args[1].args[1] != :level && __throw_log_macro_error(
            file,
            line,
            "Use `@log level=N expr` to log a message with a certain loglevel."
        )
        level = args[1].args[2]
        args = args[2:end]
    end
    level, args
end

macro log(args...)
    file, line = String(__source__.file), Int(__source__.line)
    level, msgs = log_macro_pruneargs(file, line, args)
    esc(:(
        if $(@__MODULE__).logging_enabled()
            if threadid() == 1
                with_logger($(@__MODULE__)._groebner_logger[]) do
                    @logmsg LogLevel($level) $(msgs...) _file = $file _line = $line
                end
            end
        else
            nothing
        end
    ))
end

###
# Logging memory usage

"""
    memory_logging_enabled() -> Bool

Specifies if the allocated memory information is logged in F4. If `false`, then
all memory logging is disabled, and entails no runtime overhead.

See also `@log_memory_locals`.
"""
memory_logging_enabled() = false

# Adapted from
# https://discourse.julialang.org/t/is-there-a-package-to-list-memory-consumption-of-selected-data-objects/85019/12
"""
    @log_memory_locals
    @log_memory_locals names...
    @log_memory_locals level=N names...

Logs the total allocated sizes of local variables. This does nothing when
logging is disabled in Groebner.

This may have a *significant runtime overhead*.

## Options

- If `names` argument is provided, only shows the variables present in `names`.
- If `level=N` argument is provided, then logging level `N` is used.
"""
macro log_memory_locals(names...)
    quote
        if $(@__MODULE__).logging_enabled() && $(@__MODULE__).memory_logging_enabled()
            locals = Base.@locals
            message = """
            Individual sizes (does not account for overlap):
            """
            for (name, refval) in locals
                if isempty($names) || (name in $names)
                    message *= "\t$name: $(Base.format_bytes(Base.summarysize(refval)))\n"
                end
            end
            message *=
                "Joint size: " * "$(Base.format_bytes(Base.summarysize(values(locals))))"
            @log level = -1 message
        else
            nothing
        end
    end
end
