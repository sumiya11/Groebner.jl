using Pkg;
pkg"activate --temp";
pkg"add Nemo";
using Nemo, Profile

@profile for i in 1:10000000
    Nemo.reconstruct(Nemo.ZZRingElem(rand(1:1000)), Nemo.ZZRingElem(2^30 + 3))
end
