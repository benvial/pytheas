#!/usr/bin/env bash

echo "INSTALLING SWIG"
echo '----------------'
export INSTALL_PATH=$TRAVIS_BUILD_DIR/opt
mkdir $INSTALL_PATH
cd $INSTALL_PATH
git clone https://github.com/swig/swig.git
cd swig
./autogen.sh
./configure --prefix=$INSTALL_PATH 
make
make install
cd ..
rm -rf swig
