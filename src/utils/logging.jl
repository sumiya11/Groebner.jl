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

function meta_formatter_groebner end

@noinline __throw_log_macro_error(file, line, error) =
    throw(ArgumentError("Invalid syntax for @log macro. $error"))

const _default_message_loglevel = LogLevel(0)
const _default_logger = @static if VERSION >= v"1.7.0"
    Ref{Logging.ConsoleLogger}(
        Logging.ConsoleLogger(
            Logging.Info,
            show_limited=false,
            meta_formatter=meta_formatter_groebner
        )
    )
else
    Ref{Logging.ConsoleLogger}(Logging.ConsoleLogger())
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

# Updates the global logging parameters in the Groebner module. 
function update_logger(; loglevel=nothing)
    !logging_enabled() && return nothing
    # Don't change anything if run from a worker thread.
    # Maybe throw a warning?
    threadid() != 1 && return nothing
    # if stream !== nothing
    #     prev_logger = _default_logger[]
    #     _default_logger[] = Logging.ConsoleLogger(
    #         stream,
    #         prev_logger.min_level,
    #         meta_formatter=meta_formatter_groebner
    #         # NOTE: the use of this keyword argument leads to fails on Julia
    #         # v1.6
    #         # show_limited=prev_logger.show_limited
    #     )
    # end
    if loglevel !== nothing
        new_logger = @static if VERSION >= v"1.7.0"
            Logging.ConsoleLogger(stderr, loglevel, meta_formatter=meta_formatter_groebner)
        else
            Logging.ConsoleLogger(
                stderr,
                Logging.LogLevel(loglevel);
                meta_formatter=meta_formatter_groebner
            )
        end
        _default_logger[] = new_logger
    end
    nothing
end

# Parses the arguments to the @log macro and returns a tuple of expressions that
# would evaluate to (loglevel, msg1, msg2, ...).
function pruneargs(file, line, args)
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
    level, msgs = pruneargs(file, line, args)
    esc(:(
        if $(@__MODULE__).logging_enabled()
            if threadid() == 1
                with_logger($(@__MODULE__)._default_logger[]) do
                    @logmsg LogLevel($level) $(msgs...)
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
