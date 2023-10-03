# Internal representations of monomial orderings that are actually used in the
# computation. The user interface is defined in monomials/orderings.jl.
#
# A user-defined ordering is compiled (or, either, just transformed) into a
# corresponding internal ordering, which is then used in low-level computations.

# All internal orderings are a subtype of this
abstract type AbstractInternalOrdering end

@noinline function __throw_monomial_ordering_inconsistent(msg, var_to_index, ord)
    throw(DomainError(
        ord,
        """The given monomial ordering is inconsistent with the input.
        $msg

        Variable mapping: $var_to_index
        Ordering: $ord"""
    ))
end

# Represents internal monomial orderings restricted onto a set of variables at
# some particular indices.
#
# `IsTrivial` is a parameter set either to `true` or `false`.
# If `IsTrivial` is `true`, then the indices in the ordering cover all of the
# variables (so that monomial comparison can be implemented more efficiently).
struct _Lex{IsTrivial} <: AbstractInternalOrdering
    indices::Vector{Int}
end

struct _DegLex{IsTrivial} <: AbstractInternalOrdering
    indices::Vector{Int}
end

struct _DegRevLex{IsTrivial} <: AbstractInternalOrdering
    indices::Vector{Int}
end

variable_indices(o::T) where {T <: Union{_Lex, _DegLex, _DegRevLex}} = o.indices

nontrivialization(o::_Lex) = _Lex{false}(o.indices)
nontrivialization(o::_DegLex) = _DegLex{false}(o.indices)
nontrivialization(o::_DegRevLex) = _DegRevLex{false}(o.indices)
nontrivialization(o) = o

Base.isequal(ord1::T, ord2::T) where {T <: Union{_Lex, _DegLex, _DegRevLex}} =
    variable_indices(ord1) == variable_indices(ord2)

# Maintains a list of variable indices and their weights.
struct _WeightedOrdering{T} <: AbstractInternalOrdering
    indices::Vector{Int}
    weights::Vector{T}

    function _WeightedOrdering(indices::Vector{Int}, weights::Vector{T}) where {T}
        @assert length(indices) == length(weights)
        new{T}(indices, weights)
    end
end

variable_indices(o::_WeightedOrdering) = o.indices

Base.isequal(ord1::T, ord2::T) where {T <: _WeightedOrdering} =
    ord1.indices == ord2.indices && ord1.weights == ord2.weights

# A product of two orderings.
struct _ProductOrdering{O1 <: AbstractInternalOrdering, O2 <: AbstractInternalOrdering} <:
       AbstractInternalOrdering
    # `ord1`, `ord2` may be `_ProductOrdering`s themselves
    ord1::O1
    ord2::O2

    function _ProductOrdering(
        ord1::O1,
        ord2::O2
    ) where {O1 <: AbstractInternalOrdering, O2 <: AbstractInternalOrdering}
        new{O1, O2}(ord1, ord2)
    end
end

variable_indices(o::_ProductOrdering) =
    union(variable_indices(o.ord1), variable_indices(o.ord2))

Base.isequal(ord1::T, ord2::T) where {T <: _ProductOrdering} =
    ord1.ord1 == ord2.ord1 && ord1.ord2 == ord2.ord2

struct _MatrixOrdering{T} <: AbstractInternalOrdering
    indices::Vector{Int}
    rows::Vector{Vector{T}}

    function _MatrixOrdering(indices::Vector{Int}, rows::Vector{Vector{T}}) where {T}
        new{T}(indices, rows)
    end
end

variable_indices(o::_MatrixOrdering) = o.indices

Base.isequal(ord1::T, ord2::T) where {T <: _MatrixOrdering} =
    ord1.indices == ord2.indices && ord1.rows == ord2.rows

function check_ordering_consistency(
    var_to_index::Dict{V, Int},
    ordering::Ord,
    ring_vars::Vector{T},
    ordering_vars::Vector{U};
    part_of_a_product::Bool=false
) where {V, T, U, Ord <: AbstractMonomialOrdering}
    if length(ring_vars) < length(ordering_vars)
        __throw_monomial_ordering_inconsistent(
            "There are too few variables in the ring: $ring_vars.",
            var_to_index,
            ordering
        )
    end
    if !part_of_a_product
        if !(length(ring_vars) == length(ordering_vars))
            __throw_monomial_ordering_inconsistent(
                "The number of variables in the ordering is not equal to that of the ring: $ring_vars.",
                var_to_index,
                ordering
            )
        end
    end
    ring_vars_str = map(repr, ring_vars)
    for var in ordering.variables
        if !in(repr(var), ring_vars_str)
            __throw_monomial_ordering_inconsistent(
                "The variable $var does not seem to belong to the ring: $ring_vars.",
                var_to_index,
                ordering
            )
        end
    end
    true
end

###
# Convert the interface ordering to a corresponding internal ordering

# Lex -> _Lex
# DegLex -> _DegLex
# DegRevLex -> _DegRevLex
for (Ord, InternalOrd) in ((Lex, _Lex), (DegLex, _DegLex), (DegRevLex, _DegRevLex))
    @eval begin
        function convert_to_internal_monomial_ordering(
            var_to_index::Dict{V, Int},
            ordering::$Ord{Nothing};
            part_of_a_product=false
        ) where {V}
            $InternalOrd{true}(1:length(var_to_index))
        end

        function convert_to_internal_monomial_ordering(
            var_to_index::Dict{V, Int},
            ordering::$Ord{T};
            part_of_a_product=false
        ) where {V, T}
            ring_vars = collect(keys(var_to_index))
            check_ordering_consistency(
                var_to_index,
                ordering,
                ring_vars,
                ordering.variables,
                part_of_a_product=part_of_a_product
            )
            var_to_index_str = Dict([repr(v) => i for (v, i) in var_to_index])
            indices = [var_to_index_str[repr(v)] for v in ordering.variables]
            if indices == collect(1:length(var_to_index_str))
                $InternalOrd{true}(indices)
            else
                $InternalOrd{false}(indices)
            end
        end
    end
end

# WeightedOrdering -> _WeightedOrdering
function convert_to_internal_monomial_ordering(
    var_to_index::Dict{V, Int},
    ord::WeightedOrdering{T};
    part_of_a_product=false
) where {V, T}
    ring_vars = collect(keys(var_to_index))
    if length(ring_vars) !== length(ord.weights)
        __throw_monomial_ordering_inconsistent(
            "There are too few variables in the ring: $ring_vars.",
            var_to_index,
            ord
        )
    end
    indices = collect(1:length(ring_vars))
    if any(e -> e < 0, ord.weights)
        __throw_monomial_ordering_inconsistent(
            "Negative weights are not supported, sorry.",
            var_to_index,
            ord
        )
    end
    weights = map(UInt, ord.weights)
    _WeightedOrdering(indices, weights)
end

# ProductOrdering -> _ProductOrdering
function convert_to_internal_monomial_ordering(
    var_to_index::Dict{V, Int},
    ord::ProductOrdering{T};
    part_of_a_product=false
) where {V, T}
    internal_ord1 = convert_to_internal_monomial_ordering(
        var_to_index,
        ord.ord1,
        part_of_a_product=true
    )
    internal_ord2 = convert_to_internal_monomial_ordering(
        var_to_index,
        ord.ord2,
        part_of_a_product=true
    )
    if !isempty(intersect(variable_indices(internal_ord1), variable_indices(internal_ord2)))
        @log level = -0 """
        Two blocks of the product ordering intersect by their variables.
        Block 1: $(ord.ord1)
        Block 2: $(ord.ord2)"""
    end
    _ProductOrdering(internal_ord1, internal_ord2)
end

# MatrixOrdering -> _MatrixOrdering
function convert_to_internal_monomial_ordering(
    var_to_index::Dict{V, Int},
    ordering::MatrixOrdering;
    part_of_a_product=false
) where {V}
    available_vars = collect(keys(var_to_index))
    m = length(ordering.rows)
    n = length(ordering.rows[1])
    if !(n == length(available_vars))
        __throw_monomial_ordering_inconsistent(
            "The number of columns in the matrix must be equal to the number of variables",
            var_to_index,
            ordering
        )
    end
    indices = collect(1:length(available_vars))
    _MatrixOrdering(indices, ordering.rows)
end
