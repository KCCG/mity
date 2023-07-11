# PyPI

## Dependencies

- `twine`
- `pip-compile` from `pip-tools`
- `keyring` (not available for WSL)

### keyring

Setting testpypi keyring:

```
keyring set 
```

Setting pypi keyring:
```
keyring set 
```

## Pre-Build Checklist

- Update version number in `mitylib/_version.py` and `verchew.ini`

## Build

Build scripts are available in `tools/`. Use `testpypi` version for pushing to TestPyPI and `pypi` version for pushing to PyPI.

# Docker

`docker build -f Dockerfile -t mitywgs-test-0.X.X .`