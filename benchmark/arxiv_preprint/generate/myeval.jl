# Adapted from https://discourse.julialang.org/t/expression-parser/41880/7
# code by Alan R. Rogers, Professor of Anthropology, University of Utah

import Nemo

function myeval(e::Union{Expr, Symbol, Number}, transform, map)
    try
        return _myeval(e, transform, map)
    catch ex
        println("Can't parse \"$e\"")
        rethrow(ex)
    end
end

function _myeval(s::Symbol, transform, map)
    if haskey(map, s)
        return map[s]
    else
        @info "Can not find $s in $map while parsing.."
        throw(ParseException("$s"))
    end
end

function _myeval(x::Number, transform, map)
    return transform(x)
end

# To parse an expression, convert the head to a singleton
# type, so that Julia can dispatch on that type.
function _myeval(e::Expr, transform, map)
    return _myeval(Val(e.head), e.args, transform, map)
end

# Call the function named in args[1]
function _myeval(::Val{:call}, args, transform, map)
    return _myeval(Val(args[1]), args[2:end], transform, map)
end

# Addition
function _myeval(::Val{:+}, args, transform, map)
    x = 0
    for arg in args
        x += _myeval(arg, transform, map)
    end
    return x
end

# Subtraction and negation
function _myeval(::Val{:-}, args, transform, map)
    len = length(args)
    if len == 1
        return -_myeval(args[1], transform, map)
    else
        return _myeval(args[1], transform, map) - _myeval(args[2], transform, map)
    end
end

# Multiplication
function _myeval(::Val{:*}, args, transform, map)
    x = 1
    for arg in args
        x *= _myeval(arg, transform, map)
    end
    return x
end

# Division
function _myeval(::Val{:/}, args, transform, map)
    return _myeval(args[1], transform, map) // _myeval(args[2], transform, map)
end

function _myeval(::Val{://}, args, transform, map)
    return _myeval(args[1], transform, map) // _myeval(args[2], transform, map)
end

# Exponentiation
function _myeval(::Val{:^}, args, transform, map)
    return _myeval(args[1], transform, map)^Int(numerator(_myeval(args[2], identity, map)))
end
