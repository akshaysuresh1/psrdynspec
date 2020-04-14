"""
setup.py -- setup script for installation and use of packages
"""
from setuptools import setup, find_packages

__version__ = '1.0.0'

install_requires = [
        'matplotlib>3.0;python_version>"3.7"',
        'astropy>4.0;python_version>"3.7"',
        'blimpy==1.4.2',
        'numpy',
        'scipy',
        'pandas'
]

setup(name='psrdynspec',
      version=__version__,
      description='Radio transient dynamic spectra processing tools',
      install_requires=install_requires,
      url='https://github.com/sakshay721/psrdynspec',
      author='Akshay Suresh',
      author_email='as3655@cornell.edu',
      packages=find_packages()
)
