# PyPI

## Dependencies

- `twine`
- `pip-compile` from `pip-tools`
- `keyring` (not available for WSL)

### keyring

Setting testpypi keyring:

```
keyring set https://test.pypi.org/legacy/ username
```

Setting pypi keyring:
```
keyring set https://upload.pypi.org/legacy/ username
```



## Build

Build scripts are available in `tools/`. Use `testpypi` version for pushing to TestPyPI and `pypi` version for pushing to PyPI.

### Prebuild Checklist

- Update version number in `mitylib/_version.py` and `mitylib/verchew.ini`

# Docker

## Dockerfile.dev

### Prebuild Checklist

- Add/remove any pip dependencies in the dockerfile as downloading from the test pypi does not support downloading dependencies. Remember to remove this for the production dockerfile.


```
RUN pip install pandas xlsxwriter pyfastx scipy pysam <===== CHANGE THIS
RUN pip install -i https://test.pypi.org/simple/ mitywgs-test==0.X.X
```


`docker build -f Dockerfile.dev -t mitywgs-test-0.X.X .`