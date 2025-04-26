from setuptools import setup, find_packages

with open('README.md','r') as fh:
    long_description = fh.read()

setup(
    author="George P. Tiley",
    description="A package for estimating ploidy from allele frequency distributions and other data exploration tools",
    name="estploidy",
    version="0.1.0",
    license='MIT',
    url='https://github.com/gtiley/estploidy',
    py_modules = ['estploidy_cli'],
    packages=find_packages(include=["estploidy","estploidy.*"]),
    python_requires=">=3.8",
    install_requires=[
        'pandas>=2.2.2',
        'scikit-learn>=1.5.1',
        'click>=8.1.7'
    ],
    entry_points = '''
        [console_scripts]
        estploidy=estploidy_cli:cli
    '''
)