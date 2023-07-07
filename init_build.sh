#!/bin/sh

export WORKDIR=$PWD

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

cd $SCRIPT_DIR/.. #FIXME libReco.so requires k4FWCore library to coincide with ../../k4FWCore/install/lib/libk4FWCorePlugins.so

if [ ! -d k4FWCore ]; then
  # FIXME LCG100 does not support k4FWCore and HSF spackages do not support ivy-bridge :(
  cd k4FWCore
  mkdir build install
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=../../install .. && \
  make -j `getconf _NPROCESSORS_ONLN` && \
  make install
  cd $SCRIPT_DIR/..
fi


cd $SCRIPT_DIR/..

if [ ! -d k4SimGeant4 ]; then
  cd k4SimGeant4
  mkdir build install
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=../../install -DCMAKE_CXX_STANDARD=17 ..
  make -j `getconf _NPROCESSORS_ONLN`
  make install
  cd $SCRIPT_DIR/..
fi


cd $SCRIPT_DIR/..

if [ ! -d k4Gen ]; then
  cd k4Gen
  mkdir build install
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=../../install .. && \
  make -j `getconf _NPROCESSORS_ONLN` && \
  make install
  cd $SCRIPT_DIR/..
fi

cd install

export k4FWCore_DIR=$PWD
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$k4FWCore_DIR/lib/cmake/k4FWCore
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$k4FWCore_DIR/lib/cmake/k4Gen
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$k4FWCore_DIR/lib/cmake/k4SimGeant4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$k4FWCore_DIR/lib
export PYTHONPATH=$PYTHONPATH:$k4FWCore_DIR/python
export PATH=$PATH:$k4FWCore_DIR/bin

cd $WORKDIR
