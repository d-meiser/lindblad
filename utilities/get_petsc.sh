#!/bin/sh
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.5.2.tar.gz
tar xfz petsc-3.5.2.tar.gz
cd petsc-3.5.2
PETSC_ARCH=cpucmplx-rel PETSC_DIR=`pwd` ./configure \
  --with-scalar-type=complex \
  --with-mpi=0 \
  --with-shared-libraries=0 \
  --prefix=../petsc-cpucmplx \
  --with-x=0 \
  --with-ssl=0 \
  --with-fortran-kernels=0 \
  --with-pthread=0 \
  --with-mpiuni-fortran-binding=0 \
  --with-fortran-interfaces
make PETSC_DIR=/home/dmeiser/Documents/misc/qsys/build/petsc-3.5.2 PETSC_ARCH=cpucmplx-rel all
make PETSC_DIR=/home/dmeiser/Documents/misc/qsys/build/petsc-3.5.2 PETSC_ARCH=cpucmplx-rel install
cd -

