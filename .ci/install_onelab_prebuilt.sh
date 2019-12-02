#!/bin/bash

# export ONELAB_PATH=$PWD/pytheas/tools/bin
# export VERSION="stable"
# export VERSION="dev"

ONELAB_PATH=$2
VERSION="$3"

if [ $VERSION == "stable" ]; then
  export GMSH_VERSION="4.4.1"
  export GETDP_VERSION="3.2.0"
elif [ $VERSION == "dev" ]; then
  export GMSH_VERSION="git"
  export GETDP_VERSION="git"
else
  exit
fi
#
# rm -rf $ONELAB_PATH
mkdir $ONELAB_PATH
cd $ONELAB_PATH
#
#
echo "Installing onelab $VERSION for $1 in $ONELAB_PATH"
echo "-------------------------------------------------"


case $1 in
     linux)
          export GMSH_NAME=gmsh-$GMSH_VERSION-Linux64
          export GMSH_TGZ=Linux/$GMSH_NAME.tgz
          export GETDP_NAME=getdp-$GETDP_VERSION-Linux64
          export GETDP_TGZ=Linux/$GETDP_NAME\c.tgz

          ;;
     osx)
          export GMSH_NAME=gmsh-$GMSH_VERSION-MacOSX-sdk
          export GMSH_TGZ=MacOSX/$GMSH_NAME.tgz
          export GETDP_NAME=getdp-$GETDP_VERSION-MacOSX
          export GETDP_TGZ=MacOSX/$GETDP_NAME\c.tgz
          ;;
     windows)
          echo "Not supported yet."
          ;;
esac



# gmsh

echo "installing gmsh"
echo '----------------'
wget -cq http://gmsh.info/bin/$GMSH_TGZ -O gmsh.tgz
tar -xf gmsh.tgz
rm gmsh.tgz
mv $GMSH_NAME gmsh_tmp
mv gmsh_tmp/bin/gmsh $ONELAB_PATH
if [ "$1" == "osx" ]; then
    mv gmsh_tmp/lib/*.dylib $ONELAB_PATH
fi
rm -rf gmsh_tmp


# getdp

echo "installing getdp"
echo '----------------'

wget -cq http://getdp.info/bin/$GETDP_TGZ -O getdp.tgz
tar -xf getdp.tgz
rm getdp.tgz
mv $GETDP_NAME getdp_tmp
mv getdp_tmp/bin/getdp $ONELAB_PATH
rm -rf getdp_tmp
