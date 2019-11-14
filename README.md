# PynamiX

[Documentation here](https://pages.github.sydney.edu.au/scigem/pynamix/build/html/index.html)

## Installation
Work in progress. Hopefully via `pip install pynamix` but YMMV.

## Dependencies
Work in progress - should be handled in pip install. Currently requires `python3`.

## Deploying to PyPI
Run the following to make a new distribution and upload it to PyPI
```
python3 setup.py sdist
twine upload dist/*
```

## Documentation

We use `sphinx` to manage the docs. Update documentation with:
```
cd docs
make html
```
Once these are built, you can commit and push the changes to github to have them refreshed on github pages.
