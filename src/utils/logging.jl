# This file is a part of Groebner.jl. License is GNU GPL v2.

# Logging for Groebner

# Provides the macro @log, which can be used to log a record with the given
# message to console. The macro @log wraps Logging.@logmsg from the standard
# library.

"""
    @log loglevel expr

Logs a record `expr` as a message with the given logging level.

Available options for `loglevel` are `:all`, `:debug`, `:info`, `:warn`,
`:error`.

## Examples

```jldoctest
@log :info "Hello, world!"
@log :debug "Hello, world!"
```
"""
macro log end

const _loglevels_spelled_out = (:all, :debug, :matrix, :misc, :info, :warn, :no)
const _loglevel_spell_to_int = Dict(
    :all    => -1_000_000,
    :debug  => -5,
    :matrix => -3,
    :misc   => -2,
    :info   => 0,
    :warn   => 1_000,
    :no     => 1_000_000
)
const _loglevel_default = :info

const _groebner_log_lock = Ref{ReentrantLock}(ReentrantLock())

function meta_formatter_groebner end

default_logger(io, level) = Logging.ConsoleLogger(
    io,
    level,
    show_limited=false,
    meta_formatter=meta_formatter_groebner
)

function gettime()
    tm = Libc.TmStruct(Libc.TimeVal().sec)
    (hour=tm.hour, min=tm.min, sec=tm.sec)
end

###
# Groebner logger

struct GroebnerLogger <: Logging.AbstractLogger
    logger::ConsoleLogger

    GroebnerLogger() = GroebnerLogger(stderr, Logging.LogLevel(0))
    function GroebnerLogger(io, loglevel::Logging.LogLevel)
        new(default_logger(io, loglevel))
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
end

function meta_formatter_groebner(level::LogLevel, _module, group, id, file, line)
    @nospecialize
    color = Logging.default_logcolor(level)
    time = gettime()
    h, m, s = lpad(time.hour, 2, "0"), lpad(time.min, 2, "0"), lpad(time.sec, 2, "0")
    timestr = string("[", h, ":", m, ":", s, "]")
    prefix = if level >= Logging.Warn
        "Warning"
    elseif level < Logging.Warn && level >= Logging.Info
        "Info"
    else
        "Debug"
    end
    prefix = string(timestr, " ", prefix, ":")
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
function logger_update(; loglevel=nothing)
    # Do nothing if logging is disabled
    !logging_enabled() && return nothing
    # Do nothing if nothing is being updated
    loglevel === nothing && return nothing

    new_logger = Groebner.GroebnerLogger(stderr, Logging.LogLevel(loglevel))
    lock(_groebner_log_lock[])
    try
        _groebner_logger[] = new_logger
    finally
        unlock(_groebner_log_lock[])
    end

    nothing
end

###
# The @log macro

@noinline __throw_log_macro_error(file, line, error) = throw(
    ArgumentError("""$error
                  Use as `@log loglevel message`, where `loglevel` is the logging level.
                  Supported logging levels are $_loglevels_spelled_out.""")
)

const _default_message_loglevel = Logging.Info

# Parses the arguments to the @log macro and returns a tuple of expressions that
# would evaluate to (loglevel, msg1, msg2, ...).
function log_macro_pruneargs(file, line, args)
    length(args) < 2 &&
        __throw_log_macro_error(file, line, "Argument list must contain 2 arguments.")
    !(
        args[1] isa Integer ||
        (args[1] isa QuoteNode && args[1].value in _loglevels_spelled_out)
    ) && __throw_log_macro_error(
        file,
        line,
        "Invalid logging level: $(args[1]) of type $(typeof(args[1]))"
    )
    level = args[1] isa Integer ? args[1] : _loglevel_spell_to_int[args[1].value]
    args = args[2:end]
    level, args
end

macro log(args...)
    file, line = String(__source__.file), Int(__source__.line)
    level, msgs = log_macro_pruneargs(file, line, args)
    esc(
        :(
            if $(@__MODULE__).logging_enabled()
                Groebner.Logging.with_logger($(@__MODULE__)._groebner_logger[]) do
                    $(Logging).@logmsg Groebner.Logging.LogLevel($level) $(msgs...) _file =
                        $file _line = $line
                end
            else
                nothing
            end
        )
    )
end
