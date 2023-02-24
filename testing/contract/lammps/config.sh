#!/bin/bash

set -x
set -e

source $(dirname $(readlink -f ${BASH_SOURCE}))/setup.sh

mkdir -p ${build_dir}
cd ${build_dir}

cmake \
  -DCMAKE_INSTALL_PREFIX=${install_dir} \
  -DBUILD_MPI=yes \
  -DBUILD_EXE=yes \
  -DBUILD_LIB=no \
  -DBUILD_DOC=no \
  -DLAMMPS_SIZES=smallbig \
  -DPKG_ADIOS=yes \
  ${source_dir}/cmake
