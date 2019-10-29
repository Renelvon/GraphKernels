NAME=graphkernels
PIP=pip3
PYTHON=python3
SETUP=setup.py

.PHONY: all build build_ext check install

all: build

build_ext:
	$(PYTHON) $(SETUP) build_ext

build: build_ext
	$(PYTHON) $(SETUP) build

check:
	black --check --target-version=py36 -l 80 -S --exclude=graphkernels.py setup.py $(NAME)
	check-manifest
	cpplint --recursive $(NAME)/cppkernels
	pylint setup.py $(NAME)
	pyroma -n 9 .

clean:
	git clean -xfd

dist:
	$(PYTHON) $(SETUP) sdist

distclean: clean

install: build
	$(PYTHON) $(SETUP) install --user

installcheck:
	@ cd ./Tutorial && $(PYTHON) demo_grid_search.py
	@ cd ./Tutorial && $(PYTHON) demo_mutag.py

uninstall:
	$(PIP) uninstall -y $(NAME)
