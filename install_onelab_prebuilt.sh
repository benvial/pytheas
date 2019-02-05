#!/usr/bin/env bash

export ONELAB_PATH=$HOME/bin
# export VERSION="stable"
export VERSION="dev"

if [ $VERSION == "stable" ]; then
  export GMSH_VERSION="4.0.1"
  export GETDP_VERSION="3.0.2"
elif [ $VERSION == "dev" ]; then
  export GMSH_VERSION="git"
  export GETDP_VERSION="git"
else
  exit
fi

# mkdir $ONELAB_PATH
# cd $ONELAB_PATH
# mkdir bin

echo "INSTALLING GMSH"
echo '----------------'
# gmsh
wget -c http://gmsh.info/bin/Linux/gmsh-$GMSH_VERSION-Linux64.tgz -O gmsh.tgz
tar -xvf gmsh.tgz
rm gmsh.tgz
mv gmsh-$GMSH_VERSION-Linux64/bin/gmsh $ONELAB_PATH
rm -rf gmsh

echo "INSTALLING GETDP"
echo '----------------'
# getdp
wget -c http://getdp.info/bin/Linux/getdp-$GETDP_VERSION-Linux64c.tgz -O getdp.tgz
tar -xvf getdp.tgz
rm getdp.tgz
mv getdp-$GETDP_VERSION-Linux64/bin/getdp $ONELAB_PATH
rm -rf getdp
