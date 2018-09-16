#!/usr/bin/env bash

export ONELAB_PATH=$TRAVIS_BUILD_DIR/opt
mkdir $ONELAB_PATH
cd $ONELAB_PATH
mkdir bin

echo "INSTALLING GMSH"
echo '----------------'
# gmsh
wget -c http://gmsh.info/bin/Linux/gmsh-4.0.1-Linux64.tgz -O gmsh.tgz
tar -xvf gmsh.tgz
rm gmsh.tgz
mv gmsh-3.0.6-Linux64/bin/gmsh $ONELAB_PATH/bin
rm -rf gmsh

echo "INSTALLING GETDP"
echo '----------------'
# getdp
wget -c http://getdp.info/bin/Linux/getdp-3.0.2-Linux64c.tgz -O getdp.tgz
tar -xvf getdp.tgz
rm getdp.tgz
mv getdp-2.11.3-Linux64/bin/getdp $ONELAB_PATH/bin
rm -rf getdp
