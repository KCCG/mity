#!/bin/bash

# update the requirements.txt
pip-compile 
python3 setup.py sdist
python3 setup.py bdist_wheel

version=$(grep version mitylib/_version.py | cut -f2 -d = | sed 's/[\", ]//g')

twine check dist/mitywgs-${version}*

function test {
  twine upload \
    --repository-url https://test.pypi.org/legacy/ \
    --verbose \
    --non-interactive \
    -u drmjc \
    dist/mitywgs-${version}*
}

function public {
  twine upload --verbose --non-interactive -u drmjc dist/mitywgs-${version}*
}

function docker {
  docker build --tag=latest --tag=${version} --tag=drmjc/mity:latest --tag=drmjc/mity:${version} .
  docker push drmjc/mity            # equivalent to docker push drmjc/mity:latest
  docker push drmjc/mity:${version}
  docker run drmjc/mity version
}
#test
public
wait 20 && docker
