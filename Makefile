NAME=graphkernels
PIP=pip3
PYTHON=python3
SETUP=setup.py

.PHONY: all bdist_wheel build build_ext check clean dist distclean install installcheck uninstall

all: build

bdist_wheel: build_ext
	$(PYTHON) $(SETUP) bdist_wheel

build_ext:
	$(PYTHON) $(SETUP) build_ext

build: build_ext
	$(PYTHON) $(SETUP) build

check:
	black setup.py $(NAME) Tutorial
	check-manifest
	cpplint --recursive $(NAME)/cppkernels
	pylint setup.py $(NAME) Tutorial
	pyroma -n 10 .

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
