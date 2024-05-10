# Script to run GenerateClump_Favier
# 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

# 2021 MATLAB implementation by V. Angelidakis <v.angelidakis@qub.ac>
# 2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>

from CLUMP import GenerateClump_Favier
import sys
import os
sys.path.append('../../')

inputGeom = os.path.join(os.path.dirname(__file__), 'ParticleGeometries/Ellipsoid_R_2.0_1.0_1.0.stl')
N = 10
chooseDistance = 'min'
output = 'FA_Ellipsoid_2.0_1.0_1.0.txt'
outputVTK = 'FA_Ellipsoid_2.0_1.0_1.0.vtk'
visualise = True

mesh,clump=GenerateClump_Favier(inputGeom=inputGeom, N=N, chooseDistance=chooseDistance, output=output, outputVTK=outputVTK, visualise=visualise)
