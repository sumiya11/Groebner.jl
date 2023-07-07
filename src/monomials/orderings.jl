# Monomial orderings supported by this package.

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
    Lex
    Lex(variables)

Lexicographical monomial ordering, defined as follows:

\$x_1^a_1...x_n^a_n < x_1^b_1...x_n^b_n\$

if there exists \$k ∈ 1..n\$ such that \$a_j = b_j\$ for \$j\$ in \$1..k-1\$ and
\$a_k < b_k\$.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];

# Lexicographical ordering with x > y
groebner([x*y + x, x + y^2], ordering=Lex())

# Lexicographical ordering with x < y
groebner([x*y + x, x + y^2], ordering=Lex([x, y]))
```
"""
struct Lex{T} <: AbstractMonomialOrdering
    variables::Union{Vector{T}, Nothing}
    Lex() = new{Nothing}(nothing)
    function Lex(variables::Vector{T}) where {T}
        @assert !isempty(variables)
        @assert length(unique(variables)) == length(variables)
        new{T}(variables)
    end
end

"""
    DegLex
    DegLex(variables)

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

# Degree lexicographical ordering with x < y
groebner([x*y + x, x + y^2], ordering=DegLex([x, y]))
```
"""
struct DegLex{T} <: AbstractMonomialOrdering
    variables::Union{Vector{T}, Nothing}
    DegLex() = new{Nothing}(nothing)
    function DegLex(variables::Vector{T}) where {T}
        @assert !isempty(variables)
        @assert length(unique(variables)) == length(variables)
        new{T}(variables)
    end
end

"""
    DegRevLex
    DegRevLex(variables)

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

# Degree reverse lexicographical ordering on x < y
groebner([x*y + x, x + y^2], ordering=DegRevLex())
```
"""
struct DegRevLex{T} <: AbstractMonomialOrdering
    variables::Union{Vector{T}, Nothing}
    DegRevLex() = new{Nothing}(nothing)
    function DegRevLex(variables::Vector{T}) where {T}
        @assert !isempty(variables)
        @assert length(unique(variables)) == length(variables)
        new{T}(variables)
    end
end

"""
    WeightedOrdering(weights)
    WeightedOrdering(weights, variables)

Weighted monomial ordering. Only positive weights are supported.

We use the following formal description:

TODO!

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];
ord = WeightedOrdering([3, 1])
groebner([x*y + x, x + y^2], ordering=ord)
```

Weighted ordering can be used in combination 
with some other ordering (which defaults to `Lex`)

```jldoctest
ord = WeightedOrdering([3, 1], DegRevLex())
groebner([x*y + x, x + y^2], ordering=ord)
```
"""
struct WeightedOrdering{T} <: AbstractMonomialOrdering
    weights::Vector{UInt64}
    variables::Vector{T}
    function WeightedOrdering(
        weights::Vector{T},
        variables::Union{Vector{V}, Nothing}
    ) where {T <: Integer, V}
        @assert !isempty(weights)
        @assert all(>(0), weights) "Only weights are supported."
        @assert length(weights) == length(variables)
        new{V}(weights, variables)
    end
    function WeightedOrdering(weights::Vector{T}) where {T <: Integer}
        WeightedOrdering(weights, nothing)
    end
end

# Block monomial ordering.
#
# We use the following formal description.
#
# Let <_1 and <_2 be two monomial orderings.
# <_1 acts on variables 1..k, and <_2 acts on variables k+1..n.
# The block monomial ordering corresponding to <_1 and <_2 
# is defined in the following way.
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff
#   x1^a1 ... xk^ak <_1 x1^b1 ... xk^bk
# or, 
#   x1^a1 ... xk^ak = x1^b1 ... xk^bk and
#   x{k+1}^a{k+1} ... xn^an <_2 x{k+1}^b{k+1} ... xn^bn
"""
    struct BlockOrdering

Block monomial ordering.

# Example

Block ordering with the first two variables ordered by `Lex`
and the last two variables ordered by `DegRevLex`:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, w) = QQ["x", "y", "z", "w"];
ord = BlockOrdering(1:2, Lex(), 3:4, DegRevLex())
groebner([x*y + w, y*z - w], ordering=ord)
```

Recursive creation of block orderings:

```jldoctest
using Groebner, AbstractAlgebra
R, (x1,x2,x3,x4,x5,x6) = QQ["x1","x2","x3","x4","x5","x6"];
ord_3_6 = BlockOrdering(3:4, DegRevLex(), 5:6, WeightedOrdering([0, 3]))
ord = BlockOrdering(1:2, Lex(), 3:6, ord_3_6)
groebner([x1 + x2 + x3 + x4 + x5 + x6], ordering=ord)
```
"""
struct BlockOrdering{R1, R2, O1, O2} <: AbstractMonomialOrdering
    r1::R1
    r2::R2
    ord1::O1
    ord2::O2
    function BlockOrdering(
        r1::R1,
        o1::O1,
        r2::R2,
        o2::O2
    ) where {
        R1 <: AbstractRange,
        O1 <: AbstractMonomialOrdering,
        R2 <: AbstractRange,
        O2 <: AbstractMonomialOrdering
    }
        @assert isdisjoint(r1, r2) "Only disjoint ranges in the block ordering are supported. Given: $r1 and $r2"
        @assert last(r1) + 1 == first(r2) "Only contiguous ranges in the block ordering are supported. Given: $r1 and $r2"
        new{R1, R2, O1, O2}(r1, r2, o1, o2)
    end
end

# Matrix monomial ordering.
#
# We use the following formal description.
#
# Let M be a matrix r × n with rows m1,...,mr. 
# M defines the following matrix ordering.
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff, there exists k ∈ 1..n such that
#   mj * a^T = mj * b^T for j in 1..k-1 
#   mk * a^T < mk * a^T
# where mj * a^T is the dot product of mj and a^T = (a1,...,an)^T
#
"""
    struct MatrixOrdering
    
Matrix monomial ordering.

# Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y, z, w) = QQ["x", "y", "z", "w"];
# the number of columns equal to the number of variables
ord = MatrixOrdering([
    1 0  0  2;
    0 0 -1 -2;
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
        # @assert all(row -> all(>=(0), row), rows)
        new(rows)
    end

    function MatrixOrdering(mat::Matrix{T}) where {T <: Integer}
        m, n = size(mat)
        rows = [mat[i, :] for i in 1:m]
        MatrixOrdering(rows)
    end
end
