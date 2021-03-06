#
# ifeq ($(TRAVIS_OS_NAME),windows)
# 	SHELL := cmd
# else
# 	SHELL := /bin/bash
# endif

SHELL := /bin/bash

VERSION=$(shell python3 -c "import pytheas; print(pytheas.__version__)")

ONELAB_VERSION = "dev"

default:
	@echo "\"make save\"?"

tag:
	# Make sure we're on the master branch
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	@echo "Tagging v$(VERSION)..."
	git tag v$(VERSION)
	git push --tags

pipy: package
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	twine upload dist/*

package: setup.py
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	rm -f dist/*
	rm -rf pytheas/tools/bin
	python3 setup.py sdist
	python3 setup.py bdist_wheel --universal


gh:
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	@echo "Pushing to github..."
	git add -A
	@read -p "Enter commit message: " MSG; \
	git commit -a -m "$$MSG"
	git push

publish: tag pipy

test:
	pytest ./tests -s --cov=./


clean: rmonelab rmtmp
	@find . | grep -E "(__pycache__|\.pyc|\.pyo$\)" | xargs rm -rf
	@rm -rf pytheas_pip.egg-info/ build/ dist/ tmp/
	cd docs && make clean


lstmp:
	@find . -type d -name 'tmp*'


rmtmp:
	@find . -type d -name 'tmp*' | xargs rm -rf


lsonelab:
	@find . -type f -name '*.pos' -o -name '*.pre' -o -name '*.msh' -o -name '*.res'

rmonelab:
	@find . -type f -name '*.pos' -o -name '*.pre' -o -name '*.msh' -o -name '*.res' | xargs rm -f

icons:
	cd docs/assets/icons && python gen.py

lint:
	flake8 setup.py pytheas/ tests/*.py

style:
	@echo "Styling..."
	black setup.py pytheas/ tests/*.py

onelab-linux:
	bash .ci/install_onelab_prebuilt.sh linux $(PWD)/pytheas/tools/bin $(ONELAB_VERSION)

onelab-osx:
	bash .ci/install_onelab_prebuilt.sh osx $(PWD)/pytheas/tools/bin $(ONELAB_VERSION)

pyinstall:
	bash .ci/pyinstall.sh


less:
	cd docs && make less


webdoc: less
	cd docs && make clean && make html

webdoc-noplot: less
	cd docs && make clean && make html-noplot

latexpdf:
	cd docs && make latexpdf


latexpdf-noplots:
	cd docs && make latexpdf-noplots
	cp docs/_build/latex/pytheas.pdf docs/_build/html/_downloads/pytheas.pdf

deploydoc: clean webdoc latexpdf
	git add -A
	git commit -a -m "update docs"
	git checkout gh-pages
	git merge master
	git push origin gh-pages
	git checkout master


## Show html doc in browser
showdoc:
	$(BROWSER) ./docs/_build/html/index.html

post:
	bash .ci/post.sh

save: clean style gh


## Install onelab local version (stable, linux only!)
onelab-local:
	bash .ci/install_onelab_prebuilt.sh linux $(PWD)/pytheas/tools/bin stable
