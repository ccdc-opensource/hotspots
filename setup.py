from __future__ import print_function
from setuptools import setup, find_packages
import os

print(find_packages())


def read_file(filename):
    with open(os.path.join(os.path.dirname(__file__), filename)) as file:
        return file.read()

setup(
    name="hotspots",
    author="Chris Radoux, Peter Curran, Mihaela Smilova",
    author_email="current.address@unknown.invalid",
    maintainer="David Lowe",
    maintainer_email="dlowe@ccdc.cam.ac.uk",
    license="MIT",
    version="1.0.6",
    url="https://github.com/prcurran/hotspots",
    packages=find_packages(),
    include_package_data=True,
    # scripts=['src/run_hotspot.py'],
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['numpy', # Custom install with CSD Python API
                      'csd-python-api>=2.0.0',
                      'matplotlib', # Custom install with CSD Python API
                      'scipy', # Custom install with CSD Python API
                      'scikit-learn>=0',
                      'scikit-image>=0.17.2',
                      'hdbscan>=0.8.26',
                      'pandas', # Custom install with CSD Python API
                      # 'futures >=3.1.1',
                      'tqdm==4.31.1',
                      'xmltodict==0.12.0'],
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst', '*.pkl'],
        # And include any *.mol2 files found in the 'probes' package, too:
        'hotspots': ['probes/*.mol2'],
    }
)
