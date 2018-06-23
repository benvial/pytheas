#!/usr/bin/env bash

echo "INSTALLING NLOPT"
echo '----------------'
export INSTALL_PATH=$TRAVIS_BUILD_DIR/opt
mkdir $INSTALL_PATH
cd $INSTALL_PATH
wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
tar zxvf nlopt-2.4.2.tar.gz
rm nlopt-2.4.2.tar.gz
# rm -rf nlopt
mv nlopt-2.4.2 nlopt
cd nlopt
make clean
make distclean
./configure --prefix=$INSTALL_PATH --enable-shared \
            --without-guile --without-matlab --without-octave
make
make install
cd ..
rm -rf nlopt
