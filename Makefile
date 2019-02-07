SHELL := /bin/bash

VERSION=$(shell python3 -c "import pytheas; print(pytheas.__version__)")

default:
	@echo "\"make save\"?"

tag:
	# Make sure we're on the master branch
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	@echo "Tagging v$(VERSION)..."
	git tag v$(VERSION)
	git push --tags

pipy: setup.py
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	rm -f dist/*
	rm -rf pytheas/tools/bin
	python3 setup.py sdist
	python3 setup.py bdist_wheel --universal
	twine upload dist/*


gh:
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	@echo "Pushing to github..."
	git add -A
	@read -p "Enter commit message: " MSG; \
	git commit -a -m "$$MSG"
	git push

publish: tag pipy

test:
	pytest -s --cov=./
	

clean:
	@find . | grep -E "(__pycache__|\.pyc|\.pyo$\)" | xargs rm -rf
	@rm -rf pygmsh.egg-info/ build/ dist/ tmp/

lint:
	flake8 setup.py pytheas/ tests/*.py

style:
	@echo "Styling..."
	black setup.py pytheas/ tests/*.py

onelab-linux:
	bash .ci/install_onelab_prebuilt.sh linux

onelab-osx:
	bash .ci/install_onelab_prebuilt.sh osx
	
pyinstall:
	bash .ci/pyinstall.sh


citest:
	source activate testenv ; \
	pytest -s --cov=./ ; \


post:
	bash .ci/post.sh

save: clean style gh
