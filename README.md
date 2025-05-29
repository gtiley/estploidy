# estploidy
Python package for ploidy estimation and data exploration from a VCF

## Disclaimer
This package is in development and not intended for use.

## Installation

### Conda Install
Soon

### Pip Install
The package has no dependencies outside of the Anaconda Distribution. Thus, it should be possible to install directly in the base environment by cloning the directory and using `pip install .` in the root of the repo. To create an isolated conda enironment for the install, here was the one I used for testing:

```python
conda create --name estploidy python=3.11
conda install --name estploidy scikit-learn
conda install --name estploidy pandas=2.2.2
conda install --name estploidy click
```

This estploidy environement was used to generate the *requirements.txt* and *environement.yml* files. Either could be used to install dependencies with `pip install -r requirements.txt` or `conda env create -f environment.yml`, respectively.