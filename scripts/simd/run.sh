#! /usr/bin/env -S bash -e

TAG=v1.12.0-rc1

JULIA=julia_no_vectorize
rm -rf $JULIA && mkdir $JULIA
cd $JULIA
git clone --depth 1 --branch=$TAG https://github.com/JuliaLang/julia/
cd julia
git apply ../../Julia_disable_LoopVectorizePass.patch
make -j 8
cd ../../
cp Project.toml $JULIA/
./$JULIA/julia/julia --project=$JULIA -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
./$JULIA/julia/julia -e "using InteractiveUtils; @code_native debuginfo=:none ntuple(_ -> Int32(1), 16) .+ ntuple(_ -> Int32(1), 16)"
./$JULIA/julia/julia --project=$JULIA experiment.jl

JULIA=julia_default
rm -rf $JULIA && mkdir $JULIA
cd $JULIA
git clone --depth=1 --branch=$TAG https://github.com/JuliaLang/julia/
cd julia
make -j 8
cd ../../
cp Project.toml $JULIA/
./$JULIA/julia/julia --project=$JULIA -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
./$JULIA/julia/julia -e "using InteractiveUtils; @code_native debuginfo=:none ntuple(_ -> Int32(1), 16) .+ ntuple(_ -> Int32(1), 16)"
./$JULIA/julia/julia --project=$JULIA experiment.jl
