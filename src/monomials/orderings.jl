# Monomial orderings supported by this package.
#
# The orderings in this file are exported, and exist mostly for communicating
# with the user. These orderings are converted to an internal representation at
# some stage of the computation. See monomials/internal-orderings.jl

# NOTE: only global monomial orderings are supported.

# NOTE: the functions that implement comparators for monomials are specific to
# the implementation of a monomial. Thus, monomial comparators for a specific
# monomial type are defined in the file that implements this monomial type.

# All monomial orderings are subtypes of AbstractMonomialOrdering
abstract type AbstractMonomialOrdering end

# Tells the algorithm to preserve the monomial ordering defined on the
# given input polynomials, if any.
"""
    InputOrdering()

Preserves the monomial ordering defined on the input polynomials.

This is the default value for the `ordering` keyword argument in the exported
functions.

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
*Dura Lex, sed Lex.*

We use the definition from Chapter 1, Computational Commutative Algebra 1, by
Martin Kreuzer, Lorenzo Robbiano. 

DOI: https://doi.org/10.1007/978-3-540-70628-1.

## Possible Options

- `compile`: If `true`, compiles and uses a specialized monomial comparator
  function, which may speed up the runtime (default is `false`). This option is
  experimental.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# Lexicographical ordering with x > y
groebner([x*y + x, x + y^2], ordering=Lex())

# Lexicographical ordering with y > x
groebner([x*y + x, x + y^2], ordering=Lex([y, x]))

# Lexicographical ordering with x > y
# Note that both syntax -- Lex([...]) and Lex(...) -- are allowed
groebner([x*y + x, x + y^2], ordering=Lex(x, y))
```
"""
struct Lex{T} <: AbstractMonomialOrdering
    variables::Union{Vector{T}, Nothing}
    compile::Bool

    Lex(variables...; kwargs...) = Lex(collect(variables); kwargs...)
    Lex() = new{Nothing}(nothing, false)

    function Lex(variables::Vector{T}; compile::Bool=false) where {T}
        @assert !isempty(variables)
        @assert length(unique(variables)) == length(variables) "Variables in the ordering must be unique"
        new{T}(variables, compile)
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
struct DegLex{T} <: AbstractMonomialOrdering
    variables::Union{Vector{T}, Nothing}
    compile::Bool

    DegLex(variables...; kwargs...) = DegLex(collect(variables); kwargs...)
    DegLex() = new{Nothing}(nothing, false)

    function DegLex(variables::Vector{T}) where {T}
        @assert !isempty(variables)
        @assert length(unique(variables)) == length(variables)
        new{T}(variables)
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
struct DegRevLex{T} <: AbstractMonomialOrdering
    variables::Union{Vector{T}, Nothing}
    compile::Bool

    DegRevLex(variables...; kwargs...) = DegRevLex(collect(variables); kwargs...)
    DegRevLex() = new{Nothing}(nothing, false)

    function DegRevLex(variables::Vector{T}) where {T}
        @assert !isempty(variables)
        @assert length(unique(variables)) == length(variables)
        new{T}(variables)
    end
end

"""
    WeightedOrdering(weights)

Weighted monomial ordering.

*Only positive weights are supported.*

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# x has weight 3, y has weight 1
ord = WeightedOrdering([3, 1])
groebner([x*y + x, x + y^2], ordering=ord)
```
"""
struct WeightedOrdering{U, T} <: AbstractMonomialOrdering
    weights::Vector{U}
    variables::Union{Nothing, T}

    # function WeightedOrdering(variables, weights::AbstractVector)
    #     @assert !isempty(weights)
    #     @assert all(>=(0), weights) "Only nonnegative weights are supported."
    #     # @assert length(weights) == length(variables)
    #     variables === nothing && (variables = [])
    #     new{eltype(variables)}(weights, variables)
    # end

    function WeightedOrdering(weights::Vector{T}) where {T <: Integer}
        @assert all(>=(0), weights)
        weights_unsigned = map(UInt64, weights)
        new{UInt64, Nothing}(weights_unsigned, nothing)
    end
end

"""
    ProductOrdering(ord1, ord2)

Product monomial ordering. Compares monomials by `ord1`, then breaks ties by
`ord2`.

Can also be constructed with `*`.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, w) = QQ["x", "y", "z", "w"];

# Ordering with x, y > w, z
ord = ProductOrdering(DegRevLex(x, y), DegRevLex(w, z))
groebner([x*y + w, y*z - w], ordering=ord)
```

You can also use the `*` operator:

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
ord = MatrixOrdering([
    1 0  0  2;
    0 0  1  2;
    0 1  1  1;
])
groebner([x*y + w, y*z - w], ordering=ord)
```
"""
struct MatrixOrdering <: AbstractMonomialOrdering
    rows::Vector{Vector{Int}}

    function MatrixOrdering(rows::Vector{Vector{T}}) where {T <: Integer}
        @assert all(!isempty, rows)
        @assert length(unique(map(length, rows))) == 1
        new(rows)
    end

    function MatrixOrdering(mat::Matrix{T}) where {T <: Integer}
        m, n = size(mat)
        rows = [mat[i, :] for i in 1:m]
        MatrixOrdering(rows)
    end
end
