# Script to run GenerateClump_Ferellec_McDowell
# 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

# 2021 MATLAB implementation by V. Angelidakis <v.angelidakis@qub.ac>
# 2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>

from CLUMP import GenerateClump_Ferellec_McDowell
import sys
import os
sys.path.append('../../')

inputGeom = os.path.join(os.path.dirname(__file__), 'ParticleGeometries/Torus.stl')
dmin = 0.1
rmin = 0.01
rstep = 0.01
pmax = 1.0
seed = 5
output = 'FM_Torus.txt'
outputVTK = 'FM_Torus.vtk'
visualise = True

mesh,clump=GenerateClump_Ferellec_McDowell(inputGeom=inputGeom, dmin=dmin, rmin=rmin, rstep=rstep, pmax=pmax, seed=seed, output=output, outputVTK=outputVTK, visualise=visualise)
