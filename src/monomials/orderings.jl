# This file is a part of Groebner.jl. License is GNU GPL v2.

###
# Monomial orderings.

# Only global monomial orderings are supported.

# All monomial orderings are subtypes of AbstractMonomialOrdering
abstract type AbstractMonomialOrdering end

"""
    InputOrdering()

Preserves the monomial ordering defined on the input polynomials.

This is the default value for the `ordering` keyword argument.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]

# Uses the ordering `InputOrdering`, which, in this case, 
# defaults to the lexicographical ordering with x > y
groebner([x*y + x, x + y^2])
```
"""
struct InputOrdering <: AbstractMonomialOrdering end

"""
    Lex()
    Lex(variables)
    Lex(variables...)

Lexicographical monomial ordering.

We use the definition from Chapter 1, Computational Commutative Algebra 1, by
Martin Kreuzer, Lorenzo Robbiano. 

DOI: https://doi.org/10.1007/978-3-540-70628-1.

*Dura Lex, sed Lex.*

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# Lexicographical ordering with x > y
groebner([x*y + x, x + y^2], ordering=Lex())

# Lexicographical ordering with y > x
groebner([x*y + x, x + y^2], ordering=Lex([y, x]))

# Lexicographical ordering with x > y
# Both syntax are allowed -- Lex([...]) and Lex(...)
groebner([x*y + x, x + y^2], ordering=Lex(x, y))
```
"""
struct Lex{IsSimple, T} <: AbstractMonomialOrdering
    # A note on the `IsSimple` parameter:
    #
    # We say that a monomial ordering is simple if it is one of Lex, DegLex, or
    # DegRevLex, and it has "the default" order of variables: x1 > ... > xn. In
    # this case, the `IsSimple` parameter is true. We make this distinction
    # because a simple ordering can have a specialized comparator function,
    # which may be faster than a general one.
    variables::Vector{T}

    Lex(variables...) = Lex(collect(variables))

    function Lex(variables::Vector{T}) where {T}
        !(length(unique(variables)) == length(variables)) &&
            throw(DomainError("Variables in the ordering must be unique. Got $variables"))
        issimple = isempty(variables)
        new{issimple, T}(variables)
    end
end

"""
    DegLex()
    DegLex(variables)
    DegLex(variables...)

Degree lexicographical monomial ordering.

We use the definition from Chapter 1, Computational Commutative Algebra 1, by
Martin Kreuzer, Lorenzo Robbiano. 

DOI: https://doi.org/10.1007/978-3-540-70628-1.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# Degree lexicographical ordering with x > y
groebner([x*y + x, x + y^2], ordering=DegLex())

# Degree lexicographical ordering with y > x
groebner([x*y + x, x + y^2], ordering=DegLex([y, x]))
```
"""
struct DegLex{IsSimple, T} <: AbstractMonomialOrdering
    variables::Vector{T}

    DegLex(variables...) = DegLex(collect(variables))

    function DegLex(variables::Vector{T}) where {T}
        !(length(unique(variables)) == length(variables)) &&
            throw(DomainError("Variables in the ordering must be unique. Got $variables"))
        issimple = isempty(variables)
        new{issimple, T}(variables)
    end
end

"""
    DegRevLex()
    DegRevLex(variables)
    DegRevLex(variables...)

Degree reverse lexicographical monomial ordering.

We use the definition from Chapter 1, Computational Commutative Algebra 1, by
Martin Kreuzer, Lorenzo Robbiano. 

DOI: https://doi.org/10.1007/978-3-540-70628-1.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# Degree reverse lexicographical ordering with x > y
groebner([x*y + x, x + y^2], ordering=DegRevLex())

# Degree reverse lexicographical ordering with y > x
groebner([x*y + x, x + y^2], ordering=DegRevLex(y, x))
```
"""
struct DegRevLex{IsSimple, T} <: AbstractMonomialOrdering
    variables::Vector{T}

    DegRevLex(variables...) = DegRevLex(collect(variables))

    function DegRevLex(variables::Vector{T}) where {T}
        !(length(unique(variables)) == length(variables)) &&
            throw(DomainError("Variables in the ordering must be unique. Got $variables"))
        issimple = isempty(variables)
        new{issimple, T}(variables)
    end
end

"""
    WeightedOrdering(weights)

Weighted monomial ordering.

Only positive weights are supported.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# x has weight 3, y has weight 1
ord = WeightedOrdering(x => 3, y => 1)
groebner([x*y + x, x + y^2], ordering=ord)
```
"""
struct WeightedOrdering{U, T} <: AbstractMonomialOrdering
    weights::Vector{U}
    variables::Vector{T}

    WeightedOrdering(var2weight...) = WeightedOrdering(Dict(var2weight))

    function WeightedOrdering(var2weight::AbstractDict{T, U}) where {T, U <: Integer}
        variables = Vector{T}()
        weights = Vector{U}()
        for (k, v) in var2weight
            !(v >= 0) && throw(DomainError("Weights must be positive. Got $var2weight"))
            push!(variables, k)
            push!(weights, v)
        end
        weights_unsigned = map(UInt64, weights)
        !(length(unique(variables)) == length(variables)) &&
            throw(DomainError("Variables in the ordering must be unique. Got $variables"))
        new{UInt64, T}(weights_unsigned, variables)
    end
end

"""
    ProductOrdering(ord1, ord2)

Product monomial ordering. Compares by `ord1`, breaks ties by `ord2`.

Can also be constructed with `*`.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, w) = QQ["x", "y", "z", "w"];

# Ordering with x, y > w, z
ord = ProductOrdering(DegRevLex(x, y), DegRevLex(w, z))
groebner([x*y + w, y*z - w], ordering=ord)
```

It is possible to use the `*` operator:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, t) = QQ["x", "y", "z", "t"];

ord1 = Lex(t)
ord2 = DegRevLex(x, y, z)
# t >> x, y, z
ord = ord1 * ord2
groebner([x*y*z + z, t * z - 1], ordering=ord)
```
"""
struct ProductOrdering{Ord1, Ord2} <: AbstractMonomialOrdering
    ord1::Ord1
    ord2::Ord2
    function ProductOrdering(
        ord1::Ord1,
        ord2::Ord2
    ) where {Ord1 <: AbstractMonomialOrdering, Ord2 <: AbstractMonomialOrdering}
        if !isempty(intersect(ordering_variables(ord1), ordering_variables(ord2)))
            @info """
            Two blocks of the product ordering intersect by their variables.
            Block 1: $(ord1)
            Block 2: $(ord2)"""
        end
        if ordering_is_simple(ord1) || ordering_is_simple(ord2)
            throw(DomainError("Invalid monomial ordering"))
        end
        new{Ord1, Ord2}(ord1, ord2)
    end
end

function *(
    ord1::Ord1,
    ord2::Ord2
) where {Ord1 <: AbstractMonomialOrdering, Ord2 <: AbstractMonomialOrdering}
    ProductOrdering(ord1, ord2)
end

"""
    MatrixOrdering(matrix)
    MatrixOrdering(Vector{Vector})
    
Matrix monomial ordering. 

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, w) = QQ["x", "y", "z", "w"];

# the number of columns equal to the number of variables
ord = MatrixOrdering(
    [x,y,z,w],
    [
    1 0  0  2;
    0 0  1  2;
    0 1  1  1;
    ])
groebner([x*y + w, y*z - w], ordering=ord)
```
"""
struct MatrixOrdering{T} <: AbstractMonomialOrdering
    variables::Vector{T}
    rows::Vector{Vector{Int64}}

    function MatrixOrdering(variables::Vector{V}, mat::Matrix{T}) where {V, T <: Integer}
        m, n = size(mat)
        rows = [mat[i, :] for i in 1:m]
        MatrixOrdering(variables, rows)
    end

    function MatrixOrdering(variables::Vector{V}, rows::Vector{Vector{T}}) where {V, T <: Integer}
        isempty(rows) && throw(DomainError("Invalid ordering."))
        !(length(unique(map(length, rows))) == 1) && throw(DomainError("Invalid ordering."))
        !all(row -> length(row) == length(variables), rows) &&
            throw(DomainError("Invalid ordering."))
        !(length(unique(variables)) == length(variables)) &&
            throw(DomainError("Variables in the ordering must be unique. Got $variables"))
        new{V}(variables, rows)
    end
end

###
# Ordering utilities

# Returns the (ordered) array of variables that define the ordering
ordering_variables(::InputOrdering) = []
ordering_variables(ord::Union{Lex, DegLex, DegRevLex, WeightedOrdering}) = ord.variables
ordering_variables(ord::ProductOrdering) =
    union(ordering_variables(ord.ord1), ordering_variables(ord.ord2))
ordering_variables(ord::MatrixOrdering) = ord.variables

# Equality of orderings
Base.:(==)(ord1::Lex, ord2::Lex) = ordering_variables(ord1) == ordering_variables(ord2)
Base.:(==)(ord1::DegLex, ord2::DegLex) = ordering_variables(ord1) == ordering_variables(ord2)
Base.:(==)(ord1::DegRevLex, ord2::DegRevLex) = ordering_variables(ord1) == ordering_variables(ord2)
Base.:(==)(ord1::WeightedOrdering, ord2::WeightedOrdering) =
    ordering_variables(ord1) == ordering_variables(ord2) && (ord1.weights == ord2.weights)
Base.:(==)(ord1::ProductOrdering, ord2::ProductOrdering) =
    ord1.ord1 == ord2.ord1 && ord1.ord2 == ord2.ord2
Base.:(==)(ord1::MatrixOrdering, ord2::MatrixOrdering) =
    ordering_variables(ord1) == ordering_variables(ord2) && ord1.rows == ord2.rows

ordering_is_simple(::Union{Lex{T}, DegLex{T}, DegRevLex{T}}) where {T} = T
ordering_is_simple(::AbstractMonomialOrdering) = false

ordering_make_not_simple(ord::AbstractMonomialOrdering, n::Int) = ord
ordering_make_not_simple(ord::Lex{true}, n::Int) = Lex(collect(1:n))
ordering_make_not_simple(ord::DegLex{true}, n::Int) = DegLex(collect(1:n))
ordering_make_not_simple(ord::DegRevLex{true}, n::Int) = DegRevLex(collect(1:n))

# Checking consistency against a polynomial ring
function ordering_check_consistency(nvars::Int, ord::AbstractMonomialOrdering)
    isempty(ordering_variables(ord)) && return nothing
    if nvars != length(ordering_variables(ord))
        throw(DomainError("The monomial ordering is invalid."))
    end
    if !issubset(ordering_variables(ord), collect(1:nvars))
        throw(DomainError("The monomial ordering is invalid."))
    end
    nothing
end

# Transform orderings

function map_variables(vars, varmap)
    if !(Set(vars) == Set(collect(keys(varmap))))
        # Fallback to string representation
        varmap_str = Dict(string(k) => v for (k, v) in varmap)
        vars_str = map(string, vars)
        !isempty(setdiff(vars_str, collect(keys(varmap_str)))) &&
            throw(DomainError("Invalid monomial ordering."))
        return map(v -> varmap_str[v], vars_str)
    end
    map(v -> varmap[v], vars)
end

ordering_transform(ord::InputOrdering, varmap::AbstractDict) = InputOrdering()
ordering_transform(ord::Lex, varmap::AbstractDict) =
    Lex(map_variables(ordering_variables(ord), varmap))
ordering_transform(ord::DegLex, varmap::AbstractDict) =
    DegLex(map_variables(ordering_variables(ord), varmap))
ordering_transform(ord::DegRevLex, varmap::AbstractDict) =
    DegRevLex(map_variables(ordering_variables(ord), varmap))

function ordering_transform(ord::WeightedOrdering, varmap::AbstractDict)
    WeightedOrdering(Dict(map_variables(ordering_variables(ord), varmap) .=> ord.weights))
end

function ordering_transform(ord::ProductOrdering, varmap::AbstractDict)
    ProductOrdering(ordering_transform(ord.ord1, varmap), ordering_transform(ord.ord2, varmap))
end

function ordering_transform(ord::MatrixOrdering, varmap::AbstractDict)
    MatrixOrdering(map_variables(ordering_variables(ord), varmap), ord.rows)
end

# Print orderings

Base.show(io::IO, ord::AbstractMonomialOrdering) = Base.show(io, MIME("text/plain"), ord)

Base.show(io::IO, ::MIME"text/plain", ord::InputOrdering) = print(io, "InputOrdering()")

function Base.show(
    io::IO,
    ::MIME"text/plain",
    ord::Ord
) where {Ord <: Union{Lex, DegLex, DegRevLex}}
    if ordering_is_simple(ord)
        print(io, "$(nameof(Ord))()")
    else
        print(io, "$(nameof(Ord))($(join(string.(ordering_variables(ord)), ",")))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", ord::WeightedOrdering)
    tmp = ""
    for i in 1:length(ord.variables)
        tmp *= string(ord.variables[i]) * "=>" * string(ord.weights[i])
    end
    print(io, "WeightedOrdering($(tmp)")
end

function Base.show(io::IO, ::MIME"text/plain", ord::ProductOrdering)
    print(io, ord.ord1)
    print(io, " × ")
    print(io, ord.ord2)
end

function Base.show(io::IO, ::MIME"text/plain", ord::MatrixOrdering)
    print(io, "MatrixOrdering($(join(string.(ordering_variables(ord)), ",")))\n")
    for row in ord.rows
        print(io, "  ")
        print(io, "[ " * join(string.(row), " ") * " ]\n")
    end
end
