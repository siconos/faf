#!/bin/bash -fx
#OAR --project siconos
#OAR -l /nodes=1,walltime=00:10:00
source /applis/site/nix.sh

TMPDIR=${TMPDIR:-/tmp}

mkdir $TMPDIR/$USER
cd $TMPDIR/$USER

rm -rf siconos
mkdir siconos
cd siconos

#export CC=mpicc FC=mpif90 CXX=mpic++'

cmake ~/siconos -LH -DCMAKE_INSTALL_PREFIX=$HOME/.local/ -DUSER_OPTIONS_FILE=$HOME/faf/compile/build_options_release.cmake -DWITH_OPENMP=ON  -DBLAS_DIR=$HOME/.nix-profile -DBOOST_ROOT=$HOME/.nix-profile/include -DGMP_INCLUDE_DIR=$HOME/.nix-profile/include -DGMP_LIBRARY=$HOME/.nix-profile/lib/libgmp.so -DMUMPS_LIBRARY_DIRECTORY=$HOME/.local -DFCLIB_LIBRARY_DIRECTORY=$HOME/.local
make -j8
make install
cd ..
