# PynamiX

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://github.com/scigem/PynamiX/actions/workflows/tests.yml/badge.svg)](https://github.com/scigem/PynamiX/actions/workflows/tests.yml)
[![Lint](https://github.com/scigem/PynamiX/actions/workflows/lint.yml/badge.svg)](https://github.com/scigem/PynamiX/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/scigem/PynamiX/branch/main/graph/badge.svg)](https://codecov.io/gh/scigem/PynamiX)
[![Downloads](https://pepy.tech/badge/pynamix/month)](https://pepy.tech/project/pynamix)

[Documentation here](https://scigem.github.io/PynamiX/build/html/index.html), or compile it yourself following the details below.

## Installation
You can install via `pip install pynamix`.


## Developers
Clone from github and then run:
```
pip install -e .
```

### Running Tests
Run the test suite:
```bash
python -m unittest discover -s pynamix/tests -p "test_*.py"
```

See [RUN_TESTS.md](RUN_TESTS.md) for more details.

### Pre-commit Hooks
This project uses pre-commit hooks for code quality. Install them with:
```bash
pip install pre-commit
pre-commit install
pre-commit install --hook-type pre-push
```

The hooks will run:
- **black** - code formatting (on commit)
- **flake8** - linting (on commit)  
- **tests** - test suite (on push)

## Examples
Try out the included Jupyter notebook to see how to use the package.

## Dependencies
Should be handled in pip install for you. Currently requires:
- python3
- matplotlib
- numpy
- scipy
- imageio

## Documentation

We use `sphinx` to manage the docs. Update documentation with:
```
cd docs
make html
```
Once these are built, you can commit and push the changes to github to have them refreshed on github pages. You can also view them locally.

## Roadmap

A sorted implementation list is as follows:

    1. Size measurement using the FFT technique
    2. Option to choose between FFT and wavelet transform for size measurement
    3. Wrapper for James's PIV code
    4. Wrapper for James's fake radiograph generator

## Deploying to PyPI (just a reminder for Benjy, please don't try this yourself)
Run the following to make a new distribution and upload it to PyPI. **Note**: You first need to update the version number in `pyproject.toml`.
```
python3 setup.py sdist
twine upload dist/*
```
