#!/usr/bin/env bash

export ONELAB_PATH=$PWD/pytheas/tools/bin
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

rm -rf $ONELAB_PATH
mkdir $ONELAB_PATH
cd $ONELAB_PATH





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



echo "INSTALLING GMSH"
echo '----------------'
wget -c http://gmsh.info/bin/$GMSH_TGZ -O gmsh.tgz
tar -xvf gmsh.tgz
rm gmsh.tgz
mv $GMSH_NAME gmsh_tmp
mv gmsh_tmp/bin/gmsh $ONELAB_PATH
if [ "$1" == "osx" ]; then
    mv gmsh_tmp/lib/*.dylib $ONELAB_PATH
fi
rm -rf gmsh_tmp


# getdp

echo "INSTALLING GETDP"
echo '----------------'

wget -c http://getdp.info/bin/$GETDP_TGZ -O getdp.tgz
tar -xvf getdp.tgz
rm getdp.tgz
mv $GETDP_NAME getdp_tmp
mv getdp_tmp/bin/getdp $ONELAB_PATH
rm -rf getdp_tmp
