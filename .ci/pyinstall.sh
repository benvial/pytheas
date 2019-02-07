#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  sudo apt-get update
  sudo apt-get install libglu1-mesa
fi


wget -q https://repo.continuum.io/miniconda/$MINICONDA -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
hash -r
conda config --set always_yes yes --set changeps1 no
conda config --add channels conda-forge
conda update -q --all
conda info -a
conda create -q -n testenv python=3
source activate testenv
conda install --yes -n testenv cython swig pytest pytest-cov tectonic nlopt
pip install -r requirements.txt
pip install -e .
mkdir $HOME/.matplotlib && touch $HOME/.matplotlib/matplotlibrc && echo "backend:Agg" > $HOME/.matplotlib/matplotlibrc
