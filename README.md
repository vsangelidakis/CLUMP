<p align="center"><img width=50% src="https://github.com/vsangelidakis/CLUMP/blob/master/figures/CLUMP_Logo_Extended.png"></p>
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
  <a href="#architectural-features">Architectural features</a> •
  <a href="#file-tree">File tree</a> •
  <a href="#simple-example">Simple example</a> •
  <a href="#credits">Credits</a> •
  <a href="#byos-bring-your-own-scripts">BYOS</a>
</p>

---

## What CLUMP does
CLUMP is a collection of scripts to generate multi-sphere particles of overlapping or non-overlapping spheres, to approximate target geometries. The motivation behind developing CLUMP stemmed from the need to compare different clump-generation techniques, both in terms of particle morphology and mechanical performance. To this, CLUMP offers (to date) two existing and well established clump-generation techniques and proposes a new one. The generated clumps can be exported to various formats, compatible with some of the most prominent DEM codes. Last, the surface of each created clump can be extracted as a triangulated mesh, allowing for a full characterisation of particle morphology, using tools like [SHAPE](https://github.com/vsangelidakis/SHAPE).

## Architectural features
CLUMP comprises the following modules:

- __Generate_Clump__
  - Favier et al (1999)
  - Ferellec and McDowell (2010)
  - Euclidean 3D (proposed in this code)

- __Export_clump__
  - YADE
  - LAMMPS
  - EDEM
  - PFC3D

- __Characterise_Clump__
  - Surface extraction

## File tree
- __CLUMP__
  - [LICENSE](LICENSE)
  - [README.md](README.md)
  - [README.txt](README.txt)
  - __classes__ (Definition of objects)
  - __examples__
  - __figures__
  - __functions__
  - __lib__ (External dependencies)


## Simple example
This example demonstrates different approaches to generate clumps for the same target geometry. The variables below are documented within each function.

```Matlab
addpath(genpath('functions'));	% Load in-house functions
addpath(genpath('lib'));		% Load external functions (dependencies)
addpath(genpath('classes'));	% Load object-oriented architecture

% Generate clumps using the approach of Ferellec and McDowell (2010)
[mesh, clump]=GenerateClump_Ferellec_McDowell( stlFile, dmin, rmin, rstep, pmax, seed, output );

% Generate clumps using the approach proposed in this code, involving the Euclidean transform of 3D images
[mesh, clump]=GenerateClump_Euclidean_3D( stlFile, N, rMin, div, overlap, output );
```

New users are advised to start from running the available examples in the [examples](examples) folder, to get familiarised with the syntax and functionalities of CLUMP.

## Credits
CLUMP uses several external functions available within the Matlab FEX community. We want to acknowledge the following contributions:
  - Qianqian Fang - [Iso2Mesh](https://uk.mathworks.com/matlabcentral/fileexchange/68258-iso2mesh)
  - Luigi Giaccari - [Surface Reconstruction From Scattered Points Cloud](https://www.mathworks.com/matlabcentral/fileexchange/63730-surface-reconstruction-from-scattered-points-cloud)
  - Pau Micó - [stlTools](https://uk.mathworks.com/matlabcentral/fileexchange/51200-stltools)
  - Anton Semechko - [Rigid body parameters of closed surface meshes](https://uk.mathworks.com/matlabcentral/fileexchange/48913-rigid-body-parameters-of-closed-surface-meshes)

These external dependencies are added within the source code of CLUMP, to provide an out-of-the-box implementation. The licensing terms of each external dependency can be found inside the [lib](lib/) folder.

## BYOS (Bring Your Own Scripts)!
If you enjoy using CLUMP, you are welcome to require the implementation of new clump-generation approaches and features or even better contribute and share your implementations. CLUMP was created to provide a comparison of different methods, by collecting them in one place and we share this tool hoping that members of the community will find it useful. So, feel free to expand the code, propose improvements and report issues.

<h4 align="center">2020 © Vasileios Angelidakis, Sadegh Nadimi, Masahide Otsubo, Stefano Utili. <br/> Newcastle University, UK & The University of Tokyo, Japan</a></h4>