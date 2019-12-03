#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  sudo apt-get update
  sudo apt-get install libglu1-mesa
fi


if [[ $TRAVIS_OS_NAME == 'windows' ]]; then
  export EXT=exe
else
  export EXT=sh
fi

## install requirements with conda, if not found use pip
install_reqs () {
  while read requirement; do conda install --yes $requirement || pip install -U $requirement; done < $1
}


wget -q https://repo.continuum.io/miniconda/$MINICONDA -O miniconda.$EXT;
bash miniconda.$EXT -b -p $CONDA_BASE_PATH
hash -r
conda config --set always_yes yes --set changeps1 no
conda config --add channels conda-forge
conda info -a
conda create -q -n testenv python=3
source activate testenv
install_reqs requirements.txt
install_reqs requirements_test.txt
install_reqs requirements_doc.txt
pip install -e .
conda update -n testenv -q --all
mkdir $HOME/.matplotlib && touch $HOME/.matplotlib/matplotlibrc && echo "backend:Agg" > $HOME/.matplotlib/matplotlibrc
