# Logging for Groebner
#
# Provides the macro `@log`, see below.
#
# This logging is designed to be thread-safe.

const LogLevel = Int
const __global_logging_level = Ref{Atomic{LogLevel}}(Atomic{LogLevel}(0))

function update_logging_level(loglevel::LogLevel)
    atomic_xchg!(__global_logging_level[], loglevel)
    nothing 
end

macro log(expr)
    @log(LogLevel(0), expr)
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

# TODO: do thread-safe print
function _log(dir, file, line, level, msg)
    level < __logging_level[][] && return nothing
    println("[Groebner] [$level] $msg")
    nothing
end
