# The file contains descriptions of all monomial orderings
# supported by this package.

# All monomial orderings are subtypes of AbstractMonomialOrdering
abstract type AbstractMonomialOrdering end

# Lexicographic monomial ordering.
#
# Given x1^a1 x2^a2 ... xn^an and y1^b1 y2^b2 ... yn^bn, 
# the lexicographic ordering is defined as follows:
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff, there exists k = 1..n such that
#   aj = bj for j in 1..k-1 
#   ak < bk
struct Lex <: AbstractMonomialOrdering end

# Degree lexicographic monomial ordering.
#
# Given x1^a1 x2^a2 ... xn^an and y1^b1 y2^b2 ... yn^bn, 
# the lexicographic ordering is defined as follows:
#   x1^a1 x2^a2 ... xn^an < y1^b1 y2^b2 ... yn^bn
# iff,
#   sum(a1..an) < sum(b1..bn) 
# or, sum(a1..an) = sum(b1..bn),
# and there exists k = 1..n such that
#   aj = bj for j in 1..k-1 
#   ak < bk
struct DegLex <: AbstractMonomialOrdering end

# Reverse reversed lexicographic monomial ordering.
#
struct DegRevLex <: AbstractMonomialOrdering end

# Preserve the ordering from the given input polynomials, if any.
struct InputOrdering <: AbstractMonomialOrdering end

# Weighted monomial ordering.
# (currently not exported)
struct Weighted <: AbstractMonomialOrdering
    weights::Vector{UInt64}
    Weighted(weights::Vector{UInt64}) = new(weights)
end
