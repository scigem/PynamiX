import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pynamix", # Replace with your own username
    version="0.0.6",
    description='Tools for dealing with radiographs produced in DynamiX',
    url='http://github.sydney.edu.au/scigem/pynamix',
    author='Benjy Marks, François Guillard',
    author_email='benjy.marks@sydney.edu.au',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
    ],
    python_requires='>=3.6',
    install_requires=[
          'numpy',
          'matplotlib',
          'imageio',
          'scipy',
          'ipykernel',
          'ipywidgets',
          'sphinx',
          'nbsphinx',
          'pandoc',
      ],
    include_package_data=True,
)
