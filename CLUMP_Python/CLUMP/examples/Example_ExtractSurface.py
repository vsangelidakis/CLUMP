# Script to demonstrate the ExtractSurface function
# 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

# 2021 MATLAB implementation by V. Angelidakis <v.angelidakis@qub.ac>
# 2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>

import numpy as np
from CLUMP import ExtractSurface
from CLUMP.utils.STL_ReaderWriter import write_stl

clump = np.array([[1, 0, 0, 1.1],
                  [2, 1, 0, 1.1],
                  [3, 0, 0, 1.2],
                  [1, 0, 1, 1.2]])

N_sphere = 200
N_circle = 100
visualise = True

faces, vertices = ExtractSurface(clump, N_sphere, N_circle, visualise)

write_stl('Particle_surface.stl', vertices, faces)
