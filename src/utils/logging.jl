# Logging for Groebner
#
# Provides the `@log` macro, see below.

# Verbosity of logging in the Groebner module. 
const LogLevel = Int
const _default_logging_level = LogLevel(0)
const _global_logging_level = Ref{Atomic{LogLevel}}(Atomic{LogLevel}(0))

# Updates the global logging level in the Groebner module.
function update_logging_level(loglevel::LogLevel)
    atomic_xchg!(_global_logging_level[], loglevel)
    nothing
end

"""
    @log expr
    @log level=N expr

Logs a record with `expr` as a message.
Allows to specify the logging level with `level=N`.

*Examples:*

```julia
@log "Hello, world!"
@log level=1 "Hello, world!"
```
"""
macro log end

@noinline __throw_log_macro_error(file, line, error) =
    throw(LoadError(file, line, "Invalid syntax for @log macro. $error"))
function pruneargs(file, line, args)
    (length(args) > 2 || length(args) < 1) &&
        __throw_log_macro_error(file, line, "Argument list must be non-empty")
    level, message = _default_logging_level, first(args)
    length(args) == 1 && return level, message
    @assert length(args) == 2
    if args[1] isa Expr && args[1].head == :(=)
        @assert args[1].args[1] == :level
        level = args[1].args[2]
        message = args[2]
    elseif args[2] isa Expr && args[2].head == :(=)
        @assert args[2].args[1] == :level
        level = args[2].args[2]
        message = args[1]
    else
        __throw_log_macro_error(file, line, "Supported syntax is: `@log level=N message")
    end
    level, message
end

macro log(args...)
    dir = @__DIR__
    file, line = String(__source__.file), Int(__source__.line)
    level, msg = pruneargs(file, line, args)
    esc(:(
        if $(@__MODULE__).logging_enabled()
            _log($dir, $file, $line, $level, $msg)
        else
            nothing
        end
    ))
end

function _log(dir, file, line, level, msg)
    level < _global_logging_level[][] && return nothing
    if threadid() == 1
        println("[Groebner] [$level] $msg")
    end
    nothing
end
