# Adapted from PrettyNumbers.jl
# The licence is MIT.
# https://github.com/ronisbr/PrettyNumbers.jl/blob/main/LICENSE.md

"""
    pretty_number([io::IO, ::Type{String}, ]number::Number; kwargs...)

Print the `number`.

If the first argument is an `io`, then the number is printed to it. If it is
a `String`, then a string is returned with the printed number. It it is omitted,
then it defaults to `stdout`.

# Keywords

The following keywords are available to modify the number printing. Those
arguments depends on the type of the number.

## Rational

If `number` is `Rational`, then the following keywords are available:

- `compact::Bool`: If `true`, then the rational number will be printed
    compactly, in one line like `³/₄`. Otherwise, the rational number is printed
    using multiple lines, like (**Default** = `true`):


        123
        ————
        4567

## Numbers

Otherwise, the `number` is printed using the scientific notation in the base 10.
In this case, the following keywords are available:

- `always_print_base::Bool`: If `true`, then the base is always printed even if
    the base exponent is 0. (**Default** = `false`)
- `multiplication_sign::Char`: The multiplication sign that will be used between
    the significand and the decimal base, common options are `'⋅'` and `'×'`.
    (**Default** = `'x'`)
- `significand_format::String`: The format that will be used to print the
    signifcand, as described by the function [`Printf.@printf`](@ref).
    (**Default** = `"%g"`)
- `show_base::Bool`: If `true`, then the base will be printed. Otherwise, it
    will be omitted. (**Default** = `true`)
- `show_significand::Bool`: If `true`, then the significand will be printed.
    Otherwise, it will be omitted. (**Default** = `true`)
- `new_decimal_base::Union{Nothing, Number}`: If it is a number, then the
    decimal base of the number will be converted to this number. If it is
    `nothing`, then the base is not changed. (**Default** = `nothing`)

# Examples

```julia
julia> pretty_number(1906.1896)
1.90619 · 10³

julia> pretty_number(1906.1896, significand_format = "%.10f")
1.9061896000 · 10³

julia> pretty_number(1906.1896; new_decimal_base = 4)
0.190619 · 10⁴
```
"""
function pretty_number(number::Number; kwargs...)
    return pretty_number(stdout, number; kwargs...)
end

function pretty_number(io::IO, number::Number; kwargs...)
    _pn_text(io, number; kwargs...)
    return nothing
end

function pretty_number(::Type{String}, number::Number; kwargs...)
    io = IOBuffer()
    pretty_number(io, number; kwargs...)
    str = String(take!(io))
    return str
end

# Printing function for the text backend.
function _pn_text(
    io::IO,
    number::Number;
    always_print_base::Bool=false,
    compact::Bool=true,
    significand_format::String="%g",
    show_base::Bool=true,
    show_significand::Bool=true,
    multiplication_sign::Char='x',
    new_decimal_base::Union{Nothing, Number}=nothing
)
    if number isa Rational
        number_str = _render_number_text(number; compact)
    else
        number_str = _render_number_text(
            number;
            always_print_base,
            multiplication_sign,
            new_decimal_base,
            significand_format,
            show_base,
            show_significand
        )
    end

    # Output to the IO buffer
    # ==========================================================================
    print(io, number_str)

    return nothing
end

# Definition of exponents.
const _EXPONENTS = ("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
const _SUBNUMBERS = ("₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")

"""
    _get_significand_and_base(number::Number[, new_decimal_base::Integer])

Return the significand and the base of the `number` in base 10. If the parameter
`new_decimal_base` is passed, then the output base is converted to this value.
"""
function _get_significand_and_base(number::Number)
    if abs(number) > 0
        base        = floor(log10(abs(number)))
        significand = number / 10^base
        return significand, base
    else
        return number, 0
    end
end

function _get_significand_and_base(number::Number, new_decimal_base::Integer)
    if abs(number) > 0
        significand, base = _get_significand_and_base(number)
        fact = 10^(new_decimal_base - base)
        significand /= fact
        return significand, new_decimal_base
    else
        return number, new_decimal_base
    end
end

function _render_number_text(
    number::Number;
    always_print_base::Bool=false,
    multiplication_sign::Char='x',
    significand_format::String="%g",
    show_base::Bool=true,
    show_significand::Bool=true,
    new_decimal_base::Union{Nothing, Number}=nothing
)
    # Get the significand and the base.
    if new_decimal_base !== nothing
        significand, base = _get_significand_and_base(number, new_decimal_base)
    else
        significand, base = _get_significand_and_base(number)
    end

    # Significand
    # ==========================================================================

    significand_str = ""

    if show_significand
        fmt = Printf.Format(significand_format)
        significand_str = Printf.format(fmt, significand)
    end

    # Base
    # ==========================================================================

    # If `base` is 0, then only show it if the user wants.
    base_str = ""

    if show_base && ((base != 0) || always_print_base)
        # Create the string representation.
        aux = abs(base)
        exponent_str = ""

        while aux ≥ 1
            i = aux % 10
            exponent_str = _EXPONENTS[i + 1] * exponent_str
            aux = floor(aux / 10)
        end

        # If `exponent_str` is empty, then the base is 0.
        if always_print_base
            exponent_str = "⁰"
        end

        if base < 0
            exponent_str = "⁻" * exponent_str
        end

        if show_significand
            base_str *= " " * multiplication_sign * " "
        end

        base_str *= "10" * exponent_str
    end

    # Output
    # ==========================================================================

    number_str = significand_str * base_str

    return number_str
end

prettypercent(p) = string(@sprintf("%.2f", p * 100), "%")

function prettytime(t)
    if t < 1e3
        value, units = t, "ns"
    elseif t < 1e6
        value, units = t / 1e3, "μs"
    elseif t < 1e9
        value, units = t / 1e6, "ms"
    else
        value, units = t / 1e9, "s"
    end
    return string(@sprintf("%.3f", value), " ", units)
end
