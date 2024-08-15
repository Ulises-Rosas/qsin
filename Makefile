SHELL := /bin/bash

NAME    = `python3 setup.py --name`
VERSION = `python3 setup.py --version`

all: dist_dir upload

dist_dir:
	# python3 setup.py sdist bdist_wheel
	python3 -m build .

upload:
	twine upload dist/*
	rm -rf dist
	rm -rf build
	rm -rf $(NAME).egg-info

clean_test:
	rm ./test_data/test_sims/*.txt
	rm linear_overlapped_*.txt