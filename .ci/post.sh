#!/bin/bash


if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  set -e
  pip install doctr sphinx sphinx_bootstrap_theme sphinx_gallery Pillow codecov codacy-coverage
  cd docs
  make html
  make latex
  tectonic --version
  tectonic ./_build/latex/pytheas.tex
  mv ./_build/latex/pytheas.pdf ./_build/html/_downloads/pytheas.pdf
  cd ..
  doctr deploy . --built-docs docs/_build/html
  codecov
  export CODACY_PROJECT_TOKEN=50a1d9ea26e241f7bfc27117f783e4c0
  coverage xml
  python-codacy-coverage -r coverage.xml
fi
