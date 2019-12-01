#!/bin/bash


if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  source activate testenv
  set -e
  pip install -U sphinx_bootstrap_theme # update
  
  cd .ci
  
  source ./texlive/texlive_install.sh
  cd ..
  cd docs
  make html
  make latex
  cd ./_build/latex/
  # Texliveonfly will download packages automatically
  texliveonfly -c xelatex pytheas.tex
  cd ../..
  make latexpdf
  
  export PATH=/tmp/texlive/bin/x86_64-linux:$PATH
  tlmgr list --only-installed | grep -oP 'i \K.+?(?=:)'
  
  ## tectonic
  # tectonic --version
  # tectonic ./_build/latex/pytheas.tex
  mv ./_build/latex/pytheas.pdf ./_build/html/_downloads/pytheas.pdf
  cd ..
  doctr deploy . --built-docs docs/_build/html
  codecov
  export CODACY_PROJECT_TOKEN=50a1d9ea26e241f7bfc27117f783e4c0
  coverage xml
  python-codacy-coverage -r coverage.xml
fi
