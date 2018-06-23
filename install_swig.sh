#!/usr/bin/env bash

echo "INSTALLING SWIG"
echo '----------------'
git clone https://github.com/swig/swig.git
cd swig
./autogen.sh
./configure
make
sudo make install
cd ..
rm -rf swig
