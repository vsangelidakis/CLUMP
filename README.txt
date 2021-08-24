█ █ █ █ █                                                        █ █ █ █ █
 █ █ █ █         ██████╗██╗     ██╗   ██╗███╗   ███╗██████╗       █ █ █ █ 
█ █ █ █ █       ██╔════╝██║     ██║   ██║████╗ ████║██╔══██╗     █ █ █ █ █
 █ █ █ █        ██║     ██║     ██║   ██║██╔████╔██║██████╔╝      █ █ █ █ 
█ █ █ █ █       ██║     ██║     ██║   ██║██║╚██╔╝██║██╔═══╝      █ █ █ █ █
 █ █ █ █        ╚██████╗███████╗╚██████╔╝██║ ╚═╝ ██║██║           █ █ █ █ 
█ █ █ █ █        ╚═════╝╚══════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝          █ █ █ █ █
 █ █ █ █   Code Library for Universal Multi-spherical particles   █ █ █ █ 
█ █ █ █ █                                                        █ █ █ █ █


█ Contents
  • What CLUMP does
  • Architectural features
  • File tree
  • Simple example
  • Credits
  • BYOS (Bring Your Own Scripts)!

█ What CLUMP does
CLUMP is a collection of scripts to generate multi-sphere particles of overlapping or non-overlapping spheres, to approximate target geometries. The motivation behind developing CLUMP stemmed from the need to compare different clump-generation techniques, both in terms of particle morphology and mechanical performance. To this, CLUMP offers (to date) two existing and well established clump-generation techniques and proposes a new one. The generated clumps can be exported to various formats, compatible with some of the most prominent DEM codes. Last, the surface of each created clump can be extracted as a triangulated mesh, allowing for a full characterisation of particle morphology, using tools like SHAPE (https://github.com/vsangelidakis/SHAPE). 

█ Architectural features
CLUMP comprises the following modules:

• GenerateClump
  • Favier et al (1999)
  • Ferellec and McDowell (2010)
  • Euclidean 3D (proposed in this code)

• ExportClump
  • YADE
  • LAMMPS
  • EDEM
  • PFC3D

• CharacteriseClump
  • Surface extraction

█ File tree
• __CLUMP__
  • LICENSE
  • README.md (README file in "Markdown" format - This generates the html summary of CLUMP on github.com)
  • README.txt (README file in readable, unformatted format - This can be used to read the summary of CLUMP in a text editor)
  • __classes__ (Definition of objects)
  • __examples__
  • __figures__
  • __functions__
  • __lib__ (External dependencies)

█ Simple example
This example demonstrates different approaches to generate clumps for the same target geometry. The variables below are documented within each function.

addpath(genpath('functions'));	% Load in-house functions
addpath(genpath('lib'));		% Load external functions (dependencies)
addpath(genpath('classes'));	% Load object-oriented architecture

% Generate clumps using the approach of Ferellec and McDowell (2010)
[mesh, clump]=clumpGenerator_Ferellec_McDowell( stlFile, dmin, rmin, rstep, pmax, seed, output );

% Generate clumps using the approach proposed in this code, involving the Euclidean transform of 3D images
[mesh, clump]=clumpGenerator_Euclidean_3D( stlFile, N, rMin, div, overlap, output );

New users are advised to start from running the available examples in the [examples](examples) folder, to get familiarised with the syntax and functionalities of CLUMP.

█ Credits
CLUMP uses several external functions available within the Matlab FEX community. We want to acknowledge the work of the following contributions:
  • Qianqian Fang - [Iso2Mesh](https://uk.mathworks.com/matlabcentral/fileexchange/68258-iso2mesh)
  • Luigi Giaccari - [Surface Reconstruction From Scattered Points Cloud](https://www.mathworks.com/matlabcentral/fileexchange/63730-surface-reconstruction-from-scattered-points-cloud)
  • Pau Micó - [stlTools](https://uk.mathworks.com/matlabcentral/fileexchange/51200-stltools)
  • Anton Semechko - [Rigid body parameters of closed surface meshes](https://uk.mathworks.com/matlabcentral/fileexchange/48913-rigid-body-parameters-of-closed-surface-meshes)

These external dependencies are added within the source code of CLUMP, to provide an out-of-the-box implementation. The licensing terms of each external dependency can be found inside the [lib](lib/) folder.

█ BYOS (Bring Your Own Scripts)!
If you enjoy using CLUMP, you are welcome to require the implementation of new clump-generation approaches and features or even better contribute and share your implementations. CLUMP was created to provide a comparison of different methods, by collecting them in one place and we share this tool hoping that members of the community will find it useful. So, feel free to expand the code, propose improvements and report issues.

█ Copyright
2021 © Vasileios Angelidakis [1], Sadegh Nadimi [1], Masahide Otsubo [2], Stefano Utili [1]
[1] Newcastle University, UK
[2] The University of Tokyo, Japan