
# NOTE: these are internal representations of monomial orderings actually used
# in the computation. See also monomials/orderings.jl for the interface.

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

struct _Lex{Trivial} <: AbstractInternalOrdering
    indices::Vector{Int}
end

struct _DegLex{Trivial} <: AbstractInternalOrdering
    indices::Vector{Int}
end

struct _DegRevLex{Trivial} <: AbstractInternalOrdering
    indices::Vector{Int}
end

variable_indices(o::_Lex) = o.indices
variable_indices(o::_DegLex) = o.indices
variable_indices(o::_DegRevLex) = o.indices

struct _WeightedOrdering{T} <: AbstractInternalOrdering
    indices::Vector{Int}
    weights::Vector{T}

    function _WeightedOrdering(indices::Vector{Int}, weights::Vector{T}) where {T}
        @assert length(indices) == length(weights)
        new{T}(indices, weights)
    end
end

variable_indices(o::_WeightedOrdering) = o.indices

struct _ProductOrdering{O1 <: AbstractInternalOrdering, O2 <: AbstractInternalOrdering} <:
       AbstractInternalOrdering
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
    union(variable_indices(o.ord1), variable_indices(o.ord1))

struct _MatrixOrdering{T} <: AbstractInternalOrdering
    indices::Vector{Int}
    rows::Vector{Vector{T}}

    function _MatrixOrdering(indices::Vector{Int}, rows::Vector{Vector{T}}) where {T}
        new{T}(indices, rows)
    end
end

variable_indices(o::_MatrixOrdering) = o.indices

function check_ordering_ring_consistency(
    var_to_index,
    ordering::Ord,
    ring_vars::Vector{T},
    ordering_vars::Vector{U};
    part_of_a_product=false
) where {T, U, Ord <: AbstractMonomialOrdering}
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
                "The number of variables in the ordering is different from the ring: $ring_vars.",
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
            check_ordering_ring_consistency(
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
        @log level = 0 "There is an intersection of variables of two different blocks in the product ordering."
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
