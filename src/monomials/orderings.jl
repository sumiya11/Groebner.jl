# Monomial orderings supported by this package.
#
# The orderings in this file are exported, and exist solely for communicating
# with the user. These orderings are converted to an internal representation at
# some stage of the computation. See monomials/internal-orderings.jl

# NOTE: only global monomial orderings are supported.

# NOTE: the functions that actually implement comparators for monomials are
# specific to the implementation of a monomial. Thus, comparators are defined in
# the files that implement monomials.

# All monomial orderings are subtypes of AbstractMonomialOrdering
abstract type AbstractMonomialOrdering end

# Tell the algorithm to preserve the monomial ordering defined on the
# given input polynomials, if any.
"""
    InputOrdering

Preserves the monomial ordering defined on the input polynomials.
This is the default option.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]

# Uses the ordering `InputOrdering`, which, in this case, 
# defaults to the lexicographical ordering with x < y
groebner([x*y + x, x + y^2])
```
"""
struct InputOrdering <: AbstractMonomialOrdering end

"""
    Lex()
    Lex(variables)
    Lex(variables...)

Lexicographical monomial ordering, defined as follows:

\$x_1^a_1...x_n^a_n < x_1^b_1...x_n^b_n\$

if there exists \$k ∈ 1..n\$ such that \$a_j = b_j\$ for \$j\$ in \$1..k-1\$ and
\$a_k < b_k\$.

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
    DegLex
    DegLex(variables)
    DegLex(variables...)

Degree lexicographical monomial ordering, defined as follows:

\$x_1^a_1...x_n^a_n < x_1^b_1...x_n^b_n\$

- if \$a_1 +...+ a_n < b_1 +...+ b_n\$, or, 
- if \$a_1 +...+ a_n = b_1 +...+ b_n\$ and there exists \$k ∈ 1..n\$ such that
    \$a_j = b_j\$ for \$j\$ in \$1..k-1\$ and \$a_k < b_k\$.

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
    DegRevLex
    DegRevLex(variables)
    DegRevLex(variables...)

Degree reverse lexicographical monomial ordering, defined as follows:

\$x_1^a_1...x_n^a_n < x_1^b_1...x_n^b_n\$

- if \$a_1 +...+ a_n < b_1 +...+ b_n\$, or, 
- if \$a_1 +...+ a_n = b_1 +...+ b_n\$ and there exists \$k ∈ 1..n\$ such that
    \$a_j = b_j\$ for \$j\$ in \$k+1..n\$ and \$a_k > b_k\$.

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

*Only positive weights are supported at the moment.*

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

Product monomial ordering. Compares monomials by `ord1`, break ties by `ord2`.

Can also be constructed using `*`.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, w) = QQ["x", "y", "z", "w"];

# Ordering with x > y >> w > z
ord = ProductOrdering(DegRevLex(x, y), DegRevLex(w, z))
groebner([x*y + w, y*z - w], ordering=ord)
```

You can also use the `*` operator:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, t) = QQ["x", "y", "z", "t"];

ord1 = Lex(t)
ord2 = DegRevLex(x, y, z)
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

*The number of matrix columns must be equal to the number of variables.*

Let ``M`` be an ``r x n`` matrix with rows ``m_1, ..., m_r``. 
Then ``M`` defines the following ordering:

\$\$x^a < x^b\$\$

if ``m_i^T a < m_i^T b`` for some ``i``, and ``m_j^T a = m_j^T b`` for ``1 < j < i``.

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
