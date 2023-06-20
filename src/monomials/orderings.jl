# The file contains descriptions and constructors 
# of all monomial orderings supported by this package.

# Monomial comparator functions
# are specific to the monomial implementation that is being used.
# Thus, in the file with a particular monomial implementation,
# all comparators supported by the implementation are implemented.
# (see, for example, src/monoms/powervector.jl - a simple exponent vector)

# All monomial orderings are subtypes of AbstractMonomialOrdering
abstract type AbstractMonomialOrdering end

# Tell the algorithm to preserve the monomial ordering defined on the
# given input polynomials, if any.
"""
    InputOrdering

Preserves the monomial ordering defined on the input polynomials.
This is the default option.

# Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];
# uses the ordering InputOrdering, which,
# in this case, defaults to Lex
groebner([x*y + x, x + y^2])
```
"""
struct InputOrdering <: AbstractMonomialOrdering end

###
# Global orderings.

# Lexicographical monomial ordering.
#
# We use the following formal description.
#
# Given x1^a1 x2^a2 ... xn^an and y1^b1 y2^b2 ... yn^bn, 
# the lexicographical ordering is defined as follows:
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff, there exists k ∈ 1..n such that
#   aj = bj for j in 1..k-1 
#   ak < bk
"""
    Lex

Lexicographical monomial ordering.

# Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];
groebner([x*y + x, x + y^2], ordering=Lex())
```
"""
struct Lex <: AbstractMonomialOrdering end

# Degree lexicographical monomial ordering.
#
# We use the following formal description.
#
# Given x1^a1 x2^a2 ... xn^an and y1^b1 y2^b2 ... yn^bn, 
# the degree lexicographical ordering is defined as follows:
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff,
#   sum(a1..an) < sum(b1..bn) 
# or, sum(a1..an) = sum(b1..bn),
# and there exists k ∈ 1..n such that
#   aj = bj for j in 1..k-1 
#   ak < bk
"""
    struct DegLex

Degree lexicographical monomial ordering.

# Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];
groebner([x*y + x, x + y^2], ordering=DegLex())
```
"""
struct DegLex <: AbstractMonomialOrdering end

# Degree reversed lexicographical monomial ordering.
#
# We use the following formal description.
#
# Given x1^a1 x2^a2 ... xn^an and y1^b1 y2^b2 ... yn^bn, 
# the degree reversed lexicographical ordering is defined as follows:
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff,
#   sum(a1..an) < sum(b1..bn) 
# or, sum(a1..an) = sum(b1..bn),
# and there exists k ∈ 1..n such that
#   aj = bj for j in k+1..n 
#   ak > bk
"""
    struct DegRevLex

Degree reversed lexicographical monomial ordering.

# Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];
groebner([x*y + x, x + y^2], ordering=DegRevLex())
```
"""
struct DegRevLex <: AbstractMonomialOrdering end

###
# Possibly local orderings.

# Weighted monomial ordering.
#
# We use the following formal description.
#
# Let w = (w1, ..., wn) be a vector of integers. Let <* be a monomial ordering.
#
# Given w, <*, and x1^a1 x2^a2 ... xn^an and y1^b1 y2^b2 ... yn^bn, 
# the weighted orderings is defined as follows:
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff
#   w1a1 + ... + wnan < w1b1 + ... + wnbn
# or, w1a1 + ... + wnan = w1b1 + ... + wnbn, and
#   x1^a1 x2^a2 ... xn^an <* y1^b1 y2^b2 ... yn^bn
#
# In the implementation, <* defaults to Lex.
# A custom <* can be provided.
"""
    struct WeightedOrdering

Weighted monomial ordering.
Only *non-negative weights* are supported at the moment.

# Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"];
ord = WeightedOrdering([3, 1])
groebner([x*y + x, x + y^2], ordering=ord)
```

Weighted ordering can be used in combination 
with some other ordering *(which naturally defaults to `Lex`)*

```jldoctest
ord = WeightedOrdering([3, 1], DegRevLex())
groebner([x*y + x, x + y^2], ordering=ord)
```
"""
struct WeightedOrdering{O} <: AbstractMonomialOrdering
    weights::Vector{UInt64}
    ord::O
    function WeightedOrdering(
        weights::Vector{T},
        o::O
    ) where {T <: Integer, O <: AbstractMonomialOrdering}
        @assert all(>=(0), weights) "Only non-negative weights are supported at the moment."
        new{O}(weights, o)
    end
    function WeightedOrdering(weights::Vector{T}) where {T <: Integer}
        WeightedOrdering(weights, Lex())
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
