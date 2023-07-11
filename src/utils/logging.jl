# Logging for Groebner
#
# Provides the macro @log, which can be used as @log message to log a record
# with the given message. @log wraps the standard Logging.@logmsg.
# 
# Logging is disabled on all threads except the one with threadid() == 1.

const _default_logger = @static if VERSION >= v"1.8.0"
    Ref{Logging.ConsoleLogger}(Logging.ConsoleLogger(Logging.Info, show_limited=false))
else
    Ref{Logging.ConsoleLogger}(Logging.ConsoleLogger(Logging.Info))
end
const _default_message_loglevel = LogLevel(0)

# Updates the global logging parameters in the Groebner module. 
function update_logger(; stream=nothing, loglevel=nothing)
    # Don't change anything if run from a worker thread.
    # Maybe throw a warning?
    threadid() != 1 && return nothing
    if stream !== nothing
        prev_logger = _default_logger[]
        _default_logger[] = Logging.ConsoleLogger(
            stream,
            prev_logger.min_level
            # show_limited=prev_logger.show_limited
        )
    end
    if loglevel !== nothing
        prev_logger = _default_logger[]
        _default_logger[] = Logging.ConsoleLogger(
            prev_logger.stream,
            loglevel
            # show_limited=prev_logger.show_limited
        )
    end
    nothing
end

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

@noinline __throw_log_macro_error(file, line, error) =
    throw(LoadError(file, line, "Invalid syntax for @log macro. $error"))

# Parses the arguments to the @log macro and returns a tuple of expressions that
# would evaluate to (loglevel, msg1, msg2, ...).
function pruneargs(file, line, args)
    length(args) < 1 &&
        __throw_log_macro_error(file, line, "Argument list must be non-empty")
    level = _default_message_loglevel
    length(args) == 1 && return level, args
    if args[1] isa Expr && args[1].head == :(=)
        @assert args[1].args[1] == :level
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
                with_logger(_default_logger[]) do
                    @logmsg LogLevel($level) $(msgs...)
                end
            end
        else
            nothing
        end
    ))
end

function _log(level, msgs...)
    if threadid() == 1
        with_logger(_default_logger[]) do
            # @logmsg(LogLevel(level), msgs...)
        end
    end
end
