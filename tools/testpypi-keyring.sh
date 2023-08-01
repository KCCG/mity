#!/bin/bash

# update the requirements.txt
pip-compile 
python3 setup.py sdist
python3 setup.py bdist_wheel

version=$(grep version mitylib/_version.py | cut -f2 -d = | sed 's/[\", ]//g')

twine check dist/mitywgs-test-${version}*

function test {
  twine upload -r testpypi --non-interactive -u trentzz dist/mitywgs-test-${version}*
}

test