
# CLUMP

## Description

CLUMP is a collection of scripts to generate multi-sphere particles of overlapping or non-overlapping spheres, which approximate target
geometries. Multi-spheres (a.k.a. clumps) are popular in numerical simulations using the Discrete Element Method (DEM), but these geometric
object can find applications in a wider perspective, past numerical modelling. The motivation behind developing CLUMP stemmed from the need
to compare different clump-generation techniques, both in terms of particle morphology and mechanical performance using the DEM. To this
end, CLUMP offers two existing and well established clump-generation techniques and proposes a new one. The generated clumps can be
exported in various formats, compatible with some of the most prominent DEM codes. Last, the surface of each generated clump can be extracted
as a triangulated mesh and saved as an stl file, allowing for a full characterisation of particle morphology, using tools like SHAPE or 3D-printing of
physical particle replicas. CLUMP was initially developed in MATLAB and has now been fully translated in Python. Both implementations will
continue to be developed in parallel as the software evolves.


## Architectural Features
CLUMP performs the following functions:

- __GenerateClump__
  - Favier et al (1999)
  - Ferellec and McDowell (2010)
  - Euclidean 3D (2021) and Extended Euclidean 3D (proposed in this code)

- __ExportClump__
  - YADE
  - MercuryDPM
  - LAMMPS
  - EDEM
  - PFC3D

- __CharacteriseClump__
  - Surface extraction



## Installation

Running the following code in the terminal (for Linux and macOS) or in the Command Prompt (for Windows) will automatically set up the module along with its corresponding dependencies.

```bash 
pip install clump-python
```


## Quick Start

After the installation using ```pip install clump-python```, CLUMP can easily be imported as

```python
import CLUMP
```

This following examples demonstrate different approaches to generate clumps for the different target geometry. The variables below are documented within each function. 

```python
from CLUMP.examples import Example_Euclidean_3D
```

```python
from CLUMP.examples import Example_Euclidean_3D_Extended
```

```python
from CLUMP.examples import Example_Ferellec_McDowell
```

```python
from CLUMP.examples import Example_Favier
```

The following example demonstrates the surface extraction method.

```python
from CLUMP.examples import Example_ExtractSurface
```

## Simple Examples

Here is an example of generating a clump from an STL file using __Euclidean Distance Transform__ method.

```python
from CLUMP import GenerateClump_Euclidean_3D  
  
inputGeom = "path_of_STL.stl"  
N = 21  
rMin = 0  
div = 102  
overlap = 0.6  
output = 'path_of_clump_output.txt'  
outputVTK = 'path_of_clump_vtk.vtk'  
visualise = True  
  
GenerateClump_Euclidean_3D(inputGeom, N, rMin, div, overlap, output=output, outputVTK=outputVTK, visualise=visualise)
```

To run the __Extended Euclidean Transform__, one needs to specify the maximum sphere radius which will trigger the extended version. The example is as follows:

```python
from CLUMP import GenerateClump_Euclidean_3D  
  
inputGeom = "path_of_STL.stl"  
N = 21  
rMin = 0  
div = 102  
overlap = 0.6  
output = 'path_of_clump_output.txt'  
outputVTK = 'path_of_clump_vtk.vtk'  
visualise = True
rMax_ratio = 0.3 # Parameter to trigger the Extended Euclidean method
  
GenerateClump_Euclidean_3D(inputGeom, N, rMin, div, overlap, output=output, outputVTK=outputVTK, visualise=visualise, rMax_ratio=rMax_ratio)
```

Here is an example for the __Ferellec-McDowell__ clump generation method:

```python
from CLUMP import GenerateClump_Ferellec_McDowell

inputGeom = "path_of_STL.stl"  
dmin = 0.1  
rmin = 0.01  
rstep = 0.01  
pmax = 1.0  
seed = 5  
output = 'path_of_clump_output.txt'  
outputVTK = 'path_of_clump_vtk.vtk'  
visualise = True  
  
mesh,clump=GenerateClump_Ferellec_McDowell(inputGeom=inputGeom, dmin=dmin, rmin=rmin, rstep=rstep, pmax=pmax, seed=seed, output=output, outputVTK=outputVTK, visualise=visualise)
```

To generate a clump using the __Favier__ method, you can use the following example:

```python
from CLUMP import GenerateClump_Favier

inputGeom = "path_of_STL.stl"  
N = 10  
chooseDistance = 'min'  
output = 'path_of_clump_output.txt'  
outputVTK = 'path_of_clump_vtk.vtk'  
visualise = True  
  
mesh,clump=GenerateClump_Favier(inputGeom=inputGeom, N=N, chooseDistance=chooseDistance, output=output, outputVTK=outputVTK, visualise=visualise)
```

Consider a clump defined as follows:

```python
import numpy as np

clump = np.array([[1, 0, 0, 1.1],  
                  [2, 1, 0, 1.1],  
                  [3, 0, 0, 1.2],  
                  [1, 0, 1, 1.2]])
```

To extract the surface of the clump, you can use __ExtractSurface__ function as follows:

```python
from CLUMP import ExtractSurface

N_sphere = 200  
N_circle = 100  
visualise = True  
  
faces, vertices = ExtractSurface(clump, N_sphere, N_circle, visualise)
```

 ## BYOS (Bring Your Own Scripts)!

If you enjoy using CLUMP, you are welcome to require the implementation of new clump-generation approaches and features or even better contribute and share your implementations. CLUMP was created to provide a comparison of different methods, by collecting them in one place and we share this tool hoping that members of the community will find it useful. So, feel free to expand the code, propose improvements and report issues.

