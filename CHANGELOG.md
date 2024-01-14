## Groebner v0.6.4 Release notes

Now we handle both `Nemo.GF` and `Nemo.Native.GF` as the coefficient ring.

## Groebner v0.6.1-3 Release notes

Bug fixes.

## Groebner v0.6.0 Release notes

New features:
- Adds the possibility of parallel computation, can be activated by the `threaded` keyword argument in function `groebner`.
- Major speed-ups over the rationals thanks to the new learn & apply strategy and rational number reconstruction from Nemo.
- Marginal speed-ups over some finite fields due to the new implementations of arithmetic in Z/pZ.

Other:
- Added a dependency on Nemo.jl.

## Groebner v0.5.0 Release notes

New features
- Added keyword argument `statistics` to all functions in the interface. It allows one to collect and print some timings and statistics after the computation finishes.
- Added section `Development` on the documentation page.
- Now each new push is automatically benchmarked against the previous release to detect performance regressions.

A couple of bugfixes.

Versions of dependencies:
- Updated dependency on AbstractAlgebra.jl.
- Added a dependency on TimerOutputs.jl.

## Groebner v0.4.4 Release notes

Now Groebner.jl has a logo!

Versions of dependencies:
- Updated dependency on AbstractAlgebra.jl.

## Groebner v0.4.3 Release notes

Versions of dependencies:
- Added a dependency on ExprTools.jl.

Changes in the interface:

- Added keyword argument `homogenize` to the `groebner` function. If
  `homogenize=:yes` is specified, then input generators are homogenized
  (saturated accordingly) before the computation starts. This option is turned
  on by default when the monomial ordering is an elimination ordering.
- The normalization of the output bases have changed. Just as before, a basis is
  normalized to have its leading coefficients equal to one. In this release, the
  leading coefficient is selected based on the ordering of the Groebner basis
  instead of the ordering of the frontend polynomial implementation.

Other changes:

- Modified the internals of the F4 linear algebra.

- Introduced internal functions `performance_counters_enabled` and `@timed_block` for assessing performance.

## Groebner v0.4.2 Release notes

Fixed the bug with `groebner_apply!`, which crashed when called on a system over a finite field different from the one used in `groebner_learn`.

## Groebner v0.4.1 Release notes

Reinstate support for MultivariatePolynomials.jl version 0.4.

## Groebner v0.4.0 Release notes 

Minor bug fixes and improvements.

Versions of dependencies:
- Added depency on SIMD.jl.
- Updated version of MultivariatePolynomials.jl to 0.5.
- Updated version of AbstractAlgebra.jl to 0.31.

Changes in the interface:
- Keyword argument `rng` removed. Instead, the `seed` keyword argument is provided for setting the seed of internal rng.
- Changed keyword arguments `check`, `linalg`, `monoms`, and `loglevel`.
In particular, now it is possible to specify `monom = :sparse` to use a sparse monomial representation.
- Added functions `grobner_learn` and `groebner_apply!`.
