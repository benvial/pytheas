#!/bin/bash

# Build only on Linux
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  source activate testenv
  set -e
  if [  -d $HOME/.texlive  ]; then
    echo "cached Texlive found -- nothing to do";
    CACHED_TEX=1;
  else
    echo "cached Texlive not found";
  fi


  cd docs
  make clean
  make less

  unset CACHED_TEX

  if [ "$CACHED_TEX" ]; then
      echo ">>> Using cached Tex environment";
    else
      echo ">>> Building Tex environment";
      cd ../.ci
      source ./texlive/texlive_install.sh
      cd ../docs
      make latex
      cd ./_build/latex/
      # Texliveonfly will download missing packages automatically
      texliveonfly -c xelatex pytheas.tex
      cd ../..
    fi

    make latexpdf
    make html
    # Copy pdf doc to _downloads
    mv ./_build/latex/pytheas.pdf ./_build/html/_downloads/pytheas.pdf
    # Copy apt requirements doc to binder directory
    cp ./binder/apt.txt ./_build/html/binder

    # List Tex packages
    export PATH=/tmp/texlive/bin/x86_64-linux:$PATH
    tlmgr list --only-installed | grep -oP 'i \K.+?(?=:)'

    cd ..
    doctr deploy . --built-docs docs/_build/html
    codecov
    export CODACY_PROJECT_TOKEN=50a1d9ea26e241f7bfc27117f783e4c0
    coverage xml
    python-codacy-coverage -r coverage.xml
fi
