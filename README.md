<p align="center"><img width=50% src="https://github.com/vsangelidakis/CLUMP/blob/master/logo/CLUMP_Logo_Extended.png"></p>
<h2 align="center">Code Library for Universal Multi-sphere Particles</a></h2>
<p align="center">
    <a href="https://github.com/vsangelidakis/CLUMP/commits/master">
    <img src="https://img.shields.io/github/last-commit/vsangelidakis/CLUMP.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub last commit">
    <a href="https://github.com/vsangelidakis/CLUMP/issues">
    <img src="https://img.shields.io/github/issues-raw/vsangelidakis/CLUMP.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub issues">
    <a href="https://github.com/vsangelidakis/CLUMP/pulls">
    <img src="https://img.shields.io/github/issues-pr-raw/vsangelidakis/CLUMP.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub pull requests">
    <a href="https://opensource.org/licenses/GPL-3.0">
    <img src="https://img.shields.io/badge/license-GPL-blue.svg"
         alt="License">
    <a href="https://twitter.com/intent/tweet?text=Code Library for Universal Multi-sphere Particles: &url=https%3A%2F%2Fgithub.com%2Fvsangelidakis%2FCLUMP">
    <img src="https://img.shields.io/twitter/url/https/github.com/vsangelidakis/CLUMP.svg?style=flat-square&logo=twitter"
         alt="GitHub tweet">
</p>
<p align="center">
  <a href="#what-CLUMP-does">What CLUMP does</a> •
  <a href="#simple-example">Simple example</a> •
  <a href="#credits">Credits</a> •
  <a href="#byos-bring-your-own-scripts">BYOS</a>
  <a href="#acknowledging-clump">Acknowledging CLUMP</a>
</p>

---

## What CLUMP does
CLUMP is a collection of scripts to generate multi-sphere particles of overlapping or non-overlapping spheres, which approximate target geometries. Multi-spheres (a.k.a. clumps) are popular in numerical simulations using the Discrete Element Method (DEM), but these geometric object can find applications in a wider perspective, past numerical modelling. The motivation behind developing CLUMP stemmed from the need to compare different clump-generation techniques, both in terms of particle morphology and mechanical performance using the DEM. To this end, CLUMP offers two existing and well established clump-generation techniques and proposes a new one. The generated clumps can be exported in various formats, compatible with some of the most prominent DEM codes. Last, the surface of each generated clump can be extracted as a triangulated mesh and saved as an stl file, allowing for a full characterisation of particle morphology, using tools like [SHAPE](https://github.com/vsangelidakis/SHAPE) or 3D-printing of physical particle replicas. CLUMP was initially developed in MATLAB and has now been fully translated in Python. Both implementations will continue to be developed in parallel as the software evolves.

## Simple example
This example demonstrates different approaches to generate clumps for the same target geometry. The variables below are documented within each function.

MATLAB implementation of CLUMP:
```Matlab
addpath(genpath('functions'));	% Load in-house functions
addpath(genpath('lib'));	% Load external functions (dependencies)
addpath(genpath('classes'));	% Load object-oriented architecture

% Generate clumps using the approach of Ferellec and McDowell (2010)
[mesh, clump]=GenerateClump_Ferellec_McDowell( stlFile, dmin, rmin, rstep, pmax, seed, output );

% Generate clumps using the approach proposed in this code, involving the Euclidean transform of 3D images
[mesh, clump]=GenerateClump_Euclidean_3D( stlFile, N, rMin, div, overlap, output );
```

Python implementation of CLUMP:
```Python
# Generate clumps using the approach of Ferellec and McDowell (2010)
from CLUMP import GenerateClump_Ferellec_McDowell
mesh,clump=GenerateClump_Ferellec_McDowell(inputGeom, dmin, rmin, rstep, pmax, seed, output, outputVTK, visualise)

# Generate clumps using the approach proposed in this code, involving the Euclidean transform of 3D images
from CLUMP import GenerateClump_Euclidean_3D
mesh,clump = GenerateClump_Euclidean_3D(inputGeom, N, rMin, div, overlap, output, outputVTK, visualise)
```

New users are advised to start from running the available examples in the [examples](examples) folder, to get familiarised with the syntax and functionalities of CLUMP.

## Credits
The MATLAB and Python implementations of CLUMP use different sets of external dependencies.

- CLUMP_MATLAB uses several external functions available within the Matlab FEX community. We want to acknowledge the following contributions:
  - Qianqian Fang - [Iso2Mesh](https://uk.mathworks.com/matlabcentral/fileexchange/68258-iso2mesh)
  - Luigi Giaccari - [Surface Reconstruction From Scattered Points Cloud](https://www.mathworks.com/matlabcentral/fileexchange/63730-surface-reconstruction-from-scattered-points-cloud)
  - Pau Micó - [stlTools](https://uk.mathworks.com/matlabcentral/fileexchange/51200-stltools)
  - Anton Semechko - [Rigid body parameters of closed surface meshes](https://uk.mathworks.com/matlabcentral/fileexchange/48913-rigid-body-parameters-of-closed-surface-meshes)

These external dependencies are added within the source code of CLUMP, to provide an out-of-the-box implementation. The licensing terms of each external dependency can be found inside the [lib](CLUMP_MATLAB/lib/) folder.

- CLUMP_Python uses alternative dependencies that carry out the same operations as their MATLAB counterparts. Some of the dependencies have also been translated from MATLAB to Python.

  - [matplotlib](https://matplotlib.org/)
  - [numpy](https://numpy.org/)
  - [numpy-stl](https://numpy-stl.readthedocs.io/en/latest/)
  - [os](https://docs.python.org/3/library/os.html)
  - [pyvista](https://docs.pyvista.org/version/stable/)
  - [scipy](https://scipy.org/)
  - [time](https://docs.python.org/3/library/time.html)
  - [trimesh](https://trimesh.org/)

## BYOS (Bring Your Own Scripts)!
If you enjoy using CLUMP, you are welcome to request the implementation of new features or even better contribute and share your implementations of new or existing clump-generation techniques. CLUMP was created out of our intent to provide the DEM community with a means of easy comparison between different particle generation methods, by collecting them in one place and we share this tool hoping that members of the community will find it useful. So, feel free to expand the code, propose improvements and report issues.

## Acknowledging CLUMP
Angelidakis, V., Nadimi, S., Otsubo, M. and Utili, S., 2021. CLUMP: A Code Library to generate Universal Multi-sphere Particles. SoftwareX 15, p.100735.

[Download BibTeX entry](https://github.com/vsangelidakis/CLUMP/blob/master/CITATION.bib)

<h4 align="center">2021 © Vasileios Angelidakis, Sadegh Nadimi, Masahide Otsubo, Stefano Utili. <br/> Newcastle University, UK & The University of Tokyo, Japan</a></h4>
