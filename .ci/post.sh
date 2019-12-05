#!/bin/bash


if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  source activate testenv
  set -e
  pip install -U -r requirements_doc.txt

  #### sphinx-gallery from git
  pip install -U pip git+https://github.com/sphinx-gallery/sphinx-gallery.git
  ####

  cd .ci
  source ./texlive/texlive_install.sh
  cd ..
  cd docs
  make latex
  cd ./_build/latex/
  # Texliveonfly will download missing packages automatically
  texliveonfly -c xelatex pytheas.tex
  cd ../..


  make clean
  make latexpdf
  make html

  mv ./_build/latex/pytheas.pdf ./_build/html/_downloads/pytheas.pdf


  cp ./binder/apt.txt ./_build/html/binder


  ## tectonic
  # tectonic --version
  # tectonic ./_build/latex/pytheas.tex


  export PATH=/tmp/texlive/bin/x86_64-linux:$PATH
  tlmgr list --only-installed | grep -oP 'i \K.+?(?=:)'

  cd ..
  doctr deploy . --built-docs docs/_build/html

  codecov
  export CODACY_PROJECT_TOKEN=50a1d9ea26e241f7bfc27117f783e4c0
  coverage xml
  python-codacy-coverage -r coverage.xml
fi
