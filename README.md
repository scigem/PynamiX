# PynamiX

[Documentation here](https://pages.github.sydney.edu.au/scigem/PynamiX/build/html/index.html) (requires login, contact [Benjy](mailto:benjy.marks@sydney.edu.au) for access, or compile it yourself following the details below).

## Installation
Work in progress. Hopefully via `pip install pynamix` but YMMV. If installing from github, try cloning and then running:
```
python3 setup.py build
python3 setup.py install
```
If you make any changes to the source code, re-run those two lines to have your changes reflected in your installed package.

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
Run the following to make a new distribution and upload it to PyPI. **Note**: You first need to update the version number in `setup.py`.
```
python3 setup.py sdist
twine upload dist/*
```

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
