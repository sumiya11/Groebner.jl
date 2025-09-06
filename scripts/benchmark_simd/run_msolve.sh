#! /usr/bin/env -S bash -e

TAG=v0.9.0

# MPFR_BUILD=$PWD/mpfr_build

# rm -rf mpfr mpfr_build mpfr_check
# git clone --depth 1 https://gitlab.inria.fr/mpfr/mpfr.git
# cd mpfr
# ./autogen.sh
# ./configure --prefix=$MPFR_BUILD
# make
# make install

# gcc -o mpfr_check mpfr_check.c -lmpfr -lgmp
# ./mpfr_check

julia julia_to_msolve.jl

MSOLVE=msolve_no_vectorize
rm -rf $MSOLVE && mkdir $MSOLVE
cd $MSOLVE
git clone --depth 1 --branch=$TAG https://github.com/algebraic-solving/msolve
cd msolve
./autogen.sh
./configure
#--with-mpfr=$MPFR_BUILD
make
cd ../../
time ./$MSOLVE/msolve/msolve -g 2 -f -v 1 cyclic7_qq.txt -o cyclic7_qq_out.txt 
time ./$MSOLVE/msolve/msolve -g 2 -v 1 -f cyclic7_zp.txt -o cyclic7_zp_out.txt

MSOLVE=msolve_default
rm -rf $MSOLVE && mkdir $MSOLVE
cd $MSOLVE
git clone --depth 1 --branch=$TAG https://github.com/algebraic-solving/msolve
cd msolve
./autogen.sh
./configure
#--with-mpfr=$MPFR_BUILD
make
cd ../../
time ./$MSOLVE/msolve/msolve -g 2 -f -v 1 cyclic7_qq.txt -o cyclic7_qq_out.txt 
time ./$MSOLVE/msolve/msolve -g 2 -v 1 -f cyclic7_zp.txt -o cyclic7_zp_out.txt
