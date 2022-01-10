"""
setup.py -- setup script for installation and use of packages
"""
import os
from setuptools import setup, find_packages, Command
from setuptools.command import install

__version__ = '1.2.0'

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

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
      url='https://github.com/akshaysuresh1/psrdynspec',
      install_requires=install_requires,
      packages=find_packages(),
      license='BSD 3-Clause License',
      classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "License :: OSI Approved :: BSD 3-Clause License",
        "Topic :: Scientific/Engineering :: Astronomy"
      ],
      cmdclass={'clean': CleanCommand}
)
