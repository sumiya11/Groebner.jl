
@def title = "Groebner.jl — Benchmarks"
@def hasmath = false
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

# Benchmarks

Here we compare Groebner.jl vs Singular.jl groebner implementations over fields of *characteristic zero*. The table below lists runtimes for several standard benchmark systems in seconds

|   System    | Size  | Degree | Ours    | Singular |
| :---:       | :---: |  | :----: |  :---:   |
| cyclic-12   |  23   |  | 0.13s  | **0.03s**    |
| cyclic-13   |  25   |  | 0.43s  | **0.10s**    |
| katsura-9   |  145   |  | **0.70s**  | 2.76s    |
| katsura-10  |  274   |  | **5.53s**  | 25.15s    |
| noon-7      |  495   |  | **0.56s**  | 0.88s    |
| noon-8      |  1338  |  | **5.42s**  | 8.04s    |

The comparison is not entirely fair at the moment: Singular verifies obtained result, and we do not. (TODO)

 (TODO: add links and refs to systems)
 
 (TODO: compare with other CAS)
