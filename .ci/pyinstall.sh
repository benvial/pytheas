#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  sudo apt-get update
  sudo apt-get install libglu1-mesa
  sudo apt-get install node-less
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

# Check for existence of files to determine if cache exists
# If the dir doesn't exist, but is slated to be cached later,
# Travis unhelpfully creates it, which then causes "dir already exists"
# errors when you go to actually install the thing, so we must non-intuitively
# delete the file before re-creating it later.
if [ -d $HOME/miniconda/envs/testenv ]; then
  echo "cached miniconda environment found -- nothing to do";
  CACHED_ENV=1;
else
  echo "cached miniconda environment not found";
  rm -rf $CONDA_BASE_PATH;
fi

unset CACHED_ENV
# if we don't have a cached conda environment then build one, otherwise just activate the cached one
if [ "$CACHED_ENV" ]; then
    echo ">>> Using cached environment";
    source activate testenv
else
    echo ">>> Building python environment";
    wget -q https://repo.continuum.io/miniconda/$MINICONDA -O miniconda.$EXT;
    bash miniconda.$EXT -b -p $CONDA_BASE_PATH
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda config --add channels conda-forge
    conda info -a
    conda create -q -n testenv python=$PY
    source activate testenv
    install_reqs requirements.txt
    pip install -r requirements_doc.txt
    pip install -r requirements_test.txt
    conda update -n testenv -q --all
fi

pip install -e .

mkdir $HOME/.matplotlib && touch $HOME/.matplotlib/matplotlibrc && echo "backend:Agg" > $HOME/.matplotlib/matplotlibrc
