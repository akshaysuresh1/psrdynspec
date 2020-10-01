"""
setup.py -- setup script for installation and use of packages
"""
from setuptools import setup, find_packages

__version__ = '1.1.0'

install_requires = [
        'matplotlib>3.0;python_version>"3.7"',
        'astropy>=4.0;python_version>"3.7"',
        'blimpy==2.0.0',
        'numpy',
        'scipy',
        'pandas'
]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='psrdynspec',
      version=__version__,
      author='Akshay Suresh',
      author_email='as3655@cornell.edu',
      description='Python 3 tools for processing dynamic spectra of radio transients',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/sakshay721/psrdynspec',
      install_requires=install_requires,
      packages=find_packages(),
      license='BSD 3-Clause License',
      classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "License :: OSI Approved :: BSD 3-Clause License",
        "Topic :: Scientific/Engineering :: Astronomy"
      ]
)
