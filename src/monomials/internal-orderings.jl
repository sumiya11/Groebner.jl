
# NOTE: these are internal representations of monomial orderings actually used
# in the computation. See also monomials/orderings.jl for the interface.

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

###
# Convert the interface ordering to a corresponding internal ordering

# Lex -> _Lex
# DegLex -> _DegLex
# DegRevLex -> _DegRevLex
for (Ord, InternalOrd) in ((Lex, _Lex), (DegLex, _DegLex), (DegRevLex, _DegRevLex))
    @eval begin
        function convert_to_internal_monomial_ordering(
            var_to_index::Dict{V, Int},
            ord::$Ord{Nothing}
        ) where {V}
            $InternalOrd{true}(1:length(var_to_index))
        end

        function convert_to_internal_monomial_ordering(
            var_to_index::Dict{V, Int},
            ord::$Ord{T}
        ) where {V, T}
            available_vars = collect(keys(var_to_index))
            if length(available_vars) < length(ord.variables)
                __throw_monomial_ordering_inconsistent(
                    "There are too few variables in the ring: $available_vars.",
                    var_to_index,
                    ord
                )
            end
            for var in ord.variables
                if !in(repr(var), map(repr, available_vars))
                    __throw_monomial_ordering_inconsistent(
                        "The variable $var does not seem to belong to the ring: $available_vars.",
                        var_to_index,
                        ord
                    )
                end
            end
            var_to_index_str = Dict([repr(v) => i for (v, i) in var_to_index])
            indices = [var_to_index_str[repr(v)] for v in ord.variables]
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
    ord::WeightedOrdering{T}
) where {V, T}
    available_vars = collect(keys(var_to_index))
    if length(available_vars) !== length(ord.weights)
        __throw_monomial_ordering_inconsistent(
            "There are too few variables in the ring: $available_vars.",
            var_to_index,
            ord
        )
    end
    indices = collect(1:length(available_vars))
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
    ord::ProductOrdering{T}
) where {V, T}
    internal_ord1 = convert_to_internal_monomial_ordering(var_to_index, ord.ord1)
    internal_ord2 = convert_to_internal_monomial_ordering(var_to_index, ord.ord2)
    if !isempty(intersect(variable_indices(internal_ord1), variable_indices(internal_ord2)))
        @log level = 1000 "Variables in two different blocks of the product ordering intersect."
    end
    _ProductOrdering(internal_ord1, internal_ord2)
end

# MatrixOrdering -> _MatrixOrdering
function convert_to_internal_monomial_ordering(
    var_to_index::Dict{V, Int},
    ord::MatrixOrdering
) where {V}
    available_vars = collect(keys(var_to_index))
    @assert length(ord.rows[1]) == length(available_vars)
    indices = collect(1:length(available_vars))
    _MatrixOrdering(indices, ord.rows)
end
