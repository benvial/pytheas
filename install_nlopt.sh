#!/usr/bin/env bash

echo "INSTALLING NLOPT"
echo '----------------'
export INSTALL_PATH=$TRAVIS_BUILD_DIR/opt
mkdir $INSTALL_PATH
cd $INSTALL_PATH
export VERSION="2.5.0"
export NLOPT="nlopt-"$VERSION
export TAR="v$VERSION.tar.gz"
wget https://github.com/stevengj/nlopt/archive/$TAR
tar zxvf $TAR
rm $TAR
mv $NLOPT nlopt
cd nlopt
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH -DNLOPT_MATLAB=OFF -DNLOPT_OCTAVE=OFF -DNLOPT_GUILE=OFF ..
make
make install
cd ..
rm -rf nlopt
