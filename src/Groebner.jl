# This file is a part of Groebner.jl. License is GNU GPL v2.
module Groebner
# Groebner is a package for computing Gröbner bases. This is the main file.
#
# Groebner works over integers modulo a prime and over the rationals. At its
# heart, Groebner implements F4, multi-modular techniques, and tracing.
#
# Parts of Groebner were adapted from msolve:
# https://github.com/algebraic-solving/msolve
# msolve is distributed under GNU GPL v2+:
# https://github.com/algebraic-solving/msolve/blob/master/COPYING
#
# More precisely, the F4 implementation in Groebner adapts monomial hashtable
# implementation and routines for critical pair handling, symbolic
# preprocessing, and linear algebra from msolve.

###
# Global switches

"""
    invariants_enabled() -> Bool

Specifies if custom asserts and invariants are checked. If `false`, then all
checks are disabled, and entail no runtime overhead.

It is useful to enable this when debugging the Groebner package.

See also `@invariant` in `src/utils/invariants.jl`.
"""
invariants_enabled() = false

"""
    logging_enabled() -> Bool

Specifies if logging is enabled. If `false`, then all logging in Groebner is
disabled, and entails **(almost)** no runtime overhead.

See also `@log` in `src/utils/logging.jl`.
"""
logging_enabled() = true

"""
    performance_counters_enabled() -> Bool

If performance-tracking macro `@timeit` should be enabled in Groebner. 

When this is `false`, all performance counters in Groebner are disabled and
entail **(almost)** no runtime overhead.
"""
performance_counters_enabled() = false

###
# Imports

# Groebner does not provide a polynomial implementation of its own but relies on
# existing symbolic computation packages in Julia for communicating with the
# user. Groebner accepts as its input polynomials from the Julia packages
# AbstractAlgebra.jl, Nemo.jl (Oscar.jl), and MultivariatePolynomials.jl.
import AbstractAlgebra
import AbstractAlgebra: base_ring, elem_type

import Atomix

import Base: *
import Base.Threads
import Base.Threads: nthreads, threadid
import Base.MultiplicativeInverses: UnsignedMultiplicativeInverse

import Combinatorics

using ExprTools

import HostCPUFeatures:
    cpu_name, register_count, register_size, has_feature, pick_vector_width, fma_fast

using Logging

import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

# At the moment, used only for rational reconstruction
import Nemo

# For printing the tables with statistics nicely
import PrettyTables
using Printf

import Primes
import Primes: nextprime

import Random
import Random: AbstractRNG

import TimerOutputs

###
# Initialization

# Groebner may use multi-threading by default.
# 1. Set the environment variable GROEBNER_NO_THREADED to 1 to disable all
#    multi-threading in Groebner
# 2. If GROEBNER_NO_THREADED=0, the keyword argument `threaded` provided by some
#    of the functions in the interface can be used to turn on/off the threading.
const _threaded = Ref(true)

function __init__()
    _threaded[] = !(get(ENV, "GROEBNER_NO_THREADED", "") == "1")

    # Setup the global logger
    _groebner_log_lock[] = ReentrantLock()
    logger_update(loglevel=Logging.Info)

    # Setup performance counters
    _groebner_timer_lock[] = ReentrantLock()

    nothing
end

###
# Includes

include("utils/logging.jl")
include("utils/invariants.jl")
include("utils/timeit.jl")
include("utils/simd.jl")
include("utils/pretty_print.jl")
include("utils/plots.jl")
include("utils/context.jl")

# Test systems, such as katsura, cyclic, etc
include("utils/examples.jl")

# Monomial orderings.
include("monomials/orderings.jl")
include("monomials/internal_orderings.jl")

include("utils/keywords.jl")

# Monomial implementations
include("monomials/common.jl")
include("monomials/exponent_vector.jl")
include("monomials/packed_utils.jl")
include("monomials/packed_vector.jl")
include("monomials/sparse_vector.jl")

# Defines some type aliases
include("arithmetic/CompositeNumber.jl")
include("utils/types.jl")

# Fast arithmetic modulo a prime
include("arithmetic/Zp.jl")
# Not so fast arithmetic in the rationals
include("arithmetic/QQ.jl")

# Selecting algorithm parameters
include("groebner/parameters.jl")

# Input-output conversions for polynomials
include("input_output/input_output.jl")
include("input_output/AbstractAlgebra.jl")
include("input_output/DynamicPolynomials.jl")

#= generic f4 =#
include("f4/hashtable.jl")
include("f4/basis.jl")
include("f4/trace.jl")
include("f4/matrix.jl")

# Linear algebra backends
include("f4/linalg/linalg.jl")
include("f4/linalg/backend.jl")
include("f4/linalg/backend_changematrix.jl")
include("f4/linalg/backend_threaded.jl")
include("f4/linalg/backend_randomized.jl")
include("f4/linalg/backend_randomized_threaded.jl")
include("f4/linalg/backend_learn_apply.jl")
include("f4/linalg/backend_learn_apply_threaded.jl")
include("f4/linalg/backend_experimental.jl")

include("f4/sort.jl")
include("f4/f4.jl")
include("f4/learn_apply.jl")
include("reconstruction/crt.jl")
include("reconstruction/ratrec.jl")

#= more high level functions =#
include("groebner/lucky.jl")
include("groebner/state.jl")
include("groebner/correctness.jl")
include("groebner/groebner.jl")
include("groebner/groebner_with_change_matrix.jl")
include("groebner/learn_apply.jl")
include("groebner/isgroebner.jl")
include("groebner/normalform.jl")
include("groebner/autoreduce.jl")
include("groebner/homogenization.jl")

#= generic fglm implementation =#
# NOTE: this is currently not exported, and is only for internal use
include("fglm/linear.jl")
include("fglm/fglm.jl")
include("fglm/fglm_internal.jl")
include("fglm/kbase.jl")

# API
include("interface.jl")

using PrecompileTools
include("precompile.jl")

###
# Exports

export groebner, groebner_learn, groebner_apply!
export groebner_with_change_matrix
export isgroebner
export normalform

export kbase

export Lex,
    DegLex, DegRevLex, InputOrdering, WeightedOrdering, ProductOrdering, MatrixOrdering

@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Groebner

# 라헬
end # module Groebner
