# Logging for Groebner
#
# Provides the `@log` macro, see below.

# Verbosity of logging in the Groebner module. 
const LogLevel = Int
const _default_logging_level = LogLevel(0)
const _global_logging_level = Ref{Atomic{LogLevel}}(Atomic{LogLevel}(0))

# Updates the global logging level in the Groebner module.
function update_logging_level(loglevel::LogLevel)
    # the function can be called from different threads, so the atomics
    atomic_xchg!(__global_logging_level[], loglevel)
    nothing 
end

"""
    @log expr
    @log level expr

Logs a record with `expr` as a message.
Allows to specify the logging level with `level`.

*Examples:*

```julia
@log "Hello, world!"
@log LogLevel(1) "Hello, world!"
```
"""
macro log end

macro log(expr)
    @log(_default_logging_level, expr)
end

macro log(level, expr)
    esc(:(
        if $(@__MODULE__).logging_enabled()
            _log(1, 1, 1, $level, $expr)
        else
            nothing
        end
    ))
end

function _log(dir, file, line, level, msg)
    level < __global_logging_level[][] && return nothing
    var = gensym()
    :(
        $var = $(esc(msg));
        # 
        if threadid() == 1
            println("[Groebner] [$($level)] $($var)")
        end;
        nothing
    )
end
