#!/usr/bin/env bash

export ONELAB_PATH=$TRAVIS_BUILD_DIR/opt
mkdir $ONELAB_PATH
cd $ONELAB_PATH
mkdir bin

echo "INSTALLING GMSH"
echo '----------------'
# gmsh
wget -c http://gmsh.info/bin/Linux/gmsh-git-Linux64.tgz
tar -xvf gmsh-git-Linux64.tgz
rm gmsh-git-Linux64.tgz
mv gmsh-git-Linux64/bin/gmsh $ONELAB_PATH/bin
rm -rf gmsh-git-Linux64

echo "INSTALLING GETDP"
echo '----------------'
# getdp
wget -c http://getdp.info/bin/Linux/getdp-git-Linux64c.tgz
tar -xvf getdp-git-Linux64c.tgz
rm getdp-git-Linux64c.tgz
mv getdp-git-Linux64/bin/getdp $ONELAB_PATH/bin
rm -rf getdp-git-Linux64
