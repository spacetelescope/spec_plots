#!/bin/csh

## Build the source distribution (EGG).
python setup.py sdist

## Build the wheel.
python setup.py bdist_wheel --python-tag py27