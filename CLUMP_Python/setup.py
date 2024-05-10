from setuptools import setup
from setuptools import find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="clump-python",
    version="0.3.2",
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy',
        'numpy-stl',
        'pyvista',
        'scipy',
        'trimesh',
    ],
    package_data={
        # Include all *.stl and *.mat files found in any directory within the package
        '': ['**/*.stl', '**/*.mat'],
    },
    author="Utku Canbolat, Vasileios Angelidakis",
    author_email="utku.canbolat@fau.de",
    description="This Python library provides tools for creating and examining clumps using techniques: the Euclidean Distance Transform, Favier, and Ferellec-McDowell. It allows for the efficient generation of clumps and the extraction of their surfaces.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPL-3.0-only",
    keywords="Clump, Clump Generation, Euclidean Distance Transform, Favier, Ferellec-McDowell, Surface Extraction",
    url="https://github.com/vsangelidakis/CLUMP",
)
