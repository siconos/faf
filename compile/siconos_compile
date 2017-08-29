#!/bin/bash -fx
#OAR --project siconos
#OAR -l /nodes=1,walltime=00:10:00
source /applis/site/nix.sh

mkdir $TMPDIR/$USER
cd $TMPDIR/$USER

rm -rf siconos
mkdir siconos
cd siconos

compiler='CC=mpicc FC=mpif90 CXX=mpic++'
compiler=''
$compiler cmake ~/siconos -LH -DCMAKE_INSTALL_PREFIX=$HOME/.local/ -DUSER_OPTIONS_FILE=/home/acary/faf/compile/build_options_release.cmake -DWITH_OPENMP=ON  -DBLAS_DIR=$HOME/.nix-profile -DBOOST_ROOT=$HOME/.nix-profile/include -DGMP_INCLUDE_DIR=$HOME/.nix-profile/include -DGMP_LIBRARY=$HOME/.nix-profile/lib/libgmp.so
make -j8
make install
cd ..
