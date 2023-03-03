
### Benchmark systems

By default, all systems can be used either over integers modulo a prime
or over the rationals.

For each of these systems, equations are known.
*Occasionally, the origin and references are known.*

**Classic examples**.

Systems of type name-n are usually useful in benchmarking
starting from $n$ around 5.

Several of these systems are regular, in addition to being zero-dimensional.

1. cyclic-n. 
2. eco-n.
3. henrion-n.
4. katsura-n.
5. kinema. The rationals.
6. nbody-n-(sym). To be clarified.
7. noon-n.
8. reimer-n.
9. root-n.
In these systems, no S-polynomials are constructed.

**Examples from msolve** (https://gitlab.lip6.fr/eder/msolve-examples).

1. f633.
2. phuoc-1. Rationals.
3. sot-1. Rationals.
4. vor-1. Rationals.
5. KLYZ12, KLYZ13, KLYZ22, KLYZ23. It is not clear, 
if these are randomly generated or have a reference.

**Some examples from HomotopyContinuation.jl** (https://www.juliahomotopycontinuation.org/PolynomialTestSystems.jl/stable/).

Systems in this section make sense over rationals.

1. chandra-n.
2. fourbar.
3. rps-10.
4. ipp.
5. *and around 5 systems more...*.

**Examples from Gleb Pogudin** (private communication).
See also https://github.com/SciML/StructuralIdentifiability.jl/tree/master/examples.

A couple of these systems are special in a sense that input is a 100 mb file, and output is linear with simple coefficients.

1. SIWR.
2. SEAIJRC.
3. MAPK6.
4. Pharm.
5. MAPK5.
6. CRN.
7. Goodwin.
8. Fujita.

**Some examples from Ilia Ilmer** (https://github.com/sumiya11/Groebner.jl/issues/42).
See also https://github.com/alexeyovchinnikov/SIAN-Julia/tree/main/examples.

1. chol. 
2. chol-1-out.
3. chol-weighted. Same as chol, but with weigthed ordering.
4. nfkb.
5. nfkb-weighted. *-//-*.

Systems from Gleb and Ilia are based on IO-equations from identifiability assessments of ODE models.

**A couple of challenging steady-state ideals for CRNs**.

These systems are steady-state ideals of CRN models mainly from ODEbase, which *take some time* for Groebner basis computation in degree-reverse-lex ordering.

1. BIOMD0000000085.
2. BIOMD0000000086.
3. *a couple more entries...*.

**Fukuoka MQ Challenge.** (https://www.mqchallenge.org/)

System solving over small finite fields. 
6 types of systems, all of the systems are overdetermined. Probably, *only 10-20 of these systems are usable* for benchmarking Groebner bases (other ones are too hard).

**the end.**
