<p align="center"><img width=50% src="https://github.com/vsangelidakis/CLUMP/blob/master/figures/CLUMP_Logo_Extended.png"></p>
<h2 align="center">CLUmp generator for Multi-sphere Particles</a></h2>
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
    <a href="https://twitter.com/intent/tweet?text=CLUmp generator for Multi-sphere Particles: &url=https%3A%2F%2Fgithub.com%2Fvsangelidakis%2FCLUMP">
    <img src="https://img.shields.io/twitter/url/https/github.com/vsangelidakis/CLUMP.svg?style=flat-square&logo=twitter"
         alt="GitHub tweet">
</p>
<p align="center">
  <a href="#what-CLUMP-does">What CLUMP does</a> •
  <a href="#architectural-features">Architectural features</a> •
  <a href="#file-tree">File tree</a> •
  <a href="#simple-example">Simple example</a> •
  <a href="#credits">Credits</a> •
  <a href="#byos-bring-your-own-scripts">BYOS</a> •
</p>

---

## What CLUMP does
CLUMP is a collection of simple scripts to generate multi-sphere particles of overlapping or non-overlapping spheres, to approximate target geometries. The produced particle shapes can be exported to several formats, compatible with various DEM solvers.

## Architectural features
CLUMP started as an implementation of existing methods to generate clumps and clusters of spherical particles. Along the way, we developed our own methodologies, which are shared as part of the code. CLUMP supports the following approaches:

```Matlab
-Favier
-Ferellec_and_McDowell
-Euclidean_map
```

## File tree
- __CLUMP__
  - [LICENSE](LICENSE)
  - [README.md](README.md)
  - __classes__ (Definition of objects)
  - __examples__
  - __figures__
  - __functions__ (Some of the functions are overloaded as methods in the classes)
  - __lib__ (External dependencies)


## Simple example
This example demonstrates different approaches to generate clumps for the same target geometry.

```Matlab
addpath(genpath('functions'));	% Load in-house functions
addpath(genpath('lib'));	% Load external functions (dependencies)
addpath(genpath('classes'));	% Load object-oriented architecture

% Generate clumps using the approach of Favier
xxx

% Generate clumps using the approach of Ferellec and McDowell
xxx
```

New users are advised to start from running the available examples in the [examples](examples) folder, to get familiarised with the syntax and functionalities of CLUMP.

## Credits
CLUMP uses several external functions available within the Matlab FEX community. We want to acknowledge the work of the following contributions:
xxx

These external dependencies are added within the source code of CLUMP, to provide an out-of-the-box implementation. The licensing terms of each external dependency can be found inside the [lib](lib/) folder.

## BYOS (Bring Your Own Scripts)!
If you enjoy using CLUMP, you are welcome to require the implementation of new clump-generation approaches and features or even better contribute and share your implementations. CLUMP was created to provide a comparison of different methods, by collecting them in one place and we share this tool hoping that members of the community will find it useful. So, feel free to expand the code, propose improvements and report issues.

<h4 align="center">2020 © Vasileios Angelidakis, Sadegh Nadimi, Masahide Otsubo, Stefano Utili. Newcastle University, UK; The University of Tokyo, Japan</a></h4>