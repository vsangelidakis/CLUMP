import numpy as np
import CLUMP.utils.PatchNormals as PatchNormals
import CLUMP.utils.RigidBodyParameters as RigidBodyParameters
import CLUMP.utils.STL_ReaderWriter as STL_ReaderWriter
from CLUMP.utils.ParticlePlotter import clump_plotter_pyvista
from CLUMP.utils.VTK_Writer import clump_to_VTK
import trimesh

"""
Implementation of the clump-generation concept proposed by Ferellec and McDowell (2010) [1]
2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

2021 MATLAB implementation by V. Angelidakis <v.angelidakis@qub.ac>
2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>

[1] Ferellec, J.F. and McDowell, G.R., 2010. Granular Matter, 12(5), pp.459-467. DOI 10.1007/s10035-010-0205-8

The main concept of this methodology:

1. We import the surface mesh of a particle.
2. We calculate the normal of each vertex pointing inwards.
3. For a random vertex on the particle surface, we start creating tangent spheres with incremental radii along the 
   vertex normal, starting from 'rmin', with a step of 'rstep', until they meet the surface of the particle.
4. We select a new vertex randomly, which has a distance larger than 'dmin' from the existing spheres and do the same.
5. When a percentage 'pmax' of all the vertices is used to generate spheres, the generation procedure stops.
-  An optional 'seed' parameter is introduced, to generate reproducible
   clumps.
   
   
Influence of parameters
rmin:	(0,inf) Larger rmin will lead to a smaller number of spheres
dmin: [0,inf) Larger dmin will lead to a smaller number of spheres
pmax: (0,1]   Larger pmax will lead to a larger number of spheres

Pros: The authors of this methodology claim efficiency and preservation of flat faces (reduced artificial roughness 
        compared to other techniques).
Cons: The methodology is mesh-dependent, as spheres are generated at the vertices of the input mesh.
Warning: The authors of this methodology advise that if the initial mesh is very finely discretised, an adequately 
large rmin value should be used, to guard the process against "parasitic spheres", i.e. very small spheres which might 
result to numerical instabilities when using DEM.
"""


class Clump:
    def __init__(self):
        self.positions = np.empty((0, 3))
        self.radii = np.empty((0, 1))
        self.minRadius = None
        self.maxRadius = None
        self.numSpheres = None


def GenerateClump_Ferellec_McDowell(inputGeom: str, dmin: float, rmin: float, rstep: float, pmax: float, **kwargs):
    """ Function to generate clumps using the concept proposed by Ferellec and McDowell (2010)
    :param inputGeom: Directory of stl file, used to generate spheres
    :param dmin: Minimum allowed distance between new vertex of the surface mesh and existing spheres. If left zero, this distance is not cheched.
    :param rmin: Minimum radius of sphere to be generated. For coarse meshes, the actual minimum radius might be >rmin.
    :param rstep: Step used to increase the radius in each iteration, until the generated sphere meets another point of the particle.
    :param pmax: Percentage of vertices which will be used to generate spheres. The selection of vertices is random.
    :param kwargs: Can contain either of the optimal variables "seed", "output", "outputVTK", "visualise".
                - Seed value, used to achieve reproducible (random) results
                - File name for output of the clump in .txt form. If not assigned, a .txt output file is not created.
    :return: mesh: Trimesh object containing all relevant parameters of input polyhedron. The user can get:
                    mesh.vertices
                    mesh.faces
                    mesh.centroid
                    mesh.volume
                    mesh.moment_inertia
                    mesh.principal_inertia_components
                    mesh.principal_inertia_vectors

            clump: Clump object containing all relevant clump parameters
                    clump.positions		:	M-by-3 matrix containing the position of each generated sphere.
                    clump.radii			:	M-by-1 vector containing the radius of each generated sphere
                    clump.minRadius		:	Minimum generated sphere (might differ from rmin)
                    clump.maxRadius		:	Maximum generated sphere
                    clump.numSpheres	:	Total number of spheres

            output: txt file with centroids and radii, with format: [x,y,z,r]
    """

    ################################################################################################
    #                                   Main Body of the Function                                  #
    ################################################################################################

    np.random.seed(kwargs.get("seed"))  # get the argument and pass it to the random.seed function. Does nothing if None

    clump = Clump()  # instentiate Clump object for later use

    mesh = trimesh.load(inputGeom)

    F, P = mesh.faces, mesh.vertices

    # Build "mesh" structure
    mesh = RigidBodyParameters.RBP(F, P)  # calculate the rigid body parameters based on F and P.

    N = PatchNormals.patch_normals(F, P)

    Pmax = range(len(P))  # List of vertices indices

    Vertices = np.random.permutation(len(P))

    tol = rmin / 1000  # Tolerance so that the starting vertex is considered outside the sphere

    iCount = 0  # since I am stacking the arrays the counter param is not necessary
    for _ in Pmax:
        i = Vertices[iCount]
        r = rmin
        reachedMaxRadius = False

        x, y, z = P[i, 0:3]
        n = N[i, :]

        if iCount > 0 and dmin > 0:
            dcur = np.min(
                np.sqrt(np.square(x - clump.positions[:, 0].reshape(clump.positions.shape[0], 1))
                        + np.square(y - clump.positions[:, 1].reshape(clump.positions.shape[0], 1))
                        + np.square(z - clump.positions[:, 2].reshape(clump.positions.shape[0], 1)))
                - clump.radii)

            if dcur < dmin:
                iCount += 1
                continue

        while not reachedMaxRadius:
            sphMin = 1e15

            while sphMin > -tol:
                xC = x + r * n[0]
                yC = y + r * n[1]
                zC = z + r * n[2]

                distance = np.sqrt(np.square(P[:, 0] - xC)
                                   + np.square(P[:, 1] - yC)
                                   + np.square(P[:, 2] - zC))
                sph = np.square(distance / r) - 1.0
                sphMin = np.min(sph)

                r += rstep

            reachedMaxRadius = True
            indMin = np.argmin(sph) # index of the minimum

            pointInside = P[indMin, :]

            vAB = np.array([pointInside[0] - x, pointInside[1] - y, pointInside[2] - z])
            vAD = np.dot(vAB, n) / np.linalg.norm(n)

            AB = np.linalg.norm(vAB)
            AD = np.linalg.norm(vAD)

            radius = AB ** 2 / AD / 2

            xC = x + radius * n[0]
            yC = y + radius * n[1]
            zC = z + radius * n[2]

            clump.positions = np.vstack((clump.positions, np.array([xC, yC, zC]).reshape((1, 3))))
            clump.radii = np.vstack((clump.radii, radius))

        # Check whether the maximum percentage of vertices has been used
        pcur = clump.radii.shape[0] / P.shape[0]  # Current percentage of vertices used
        if pcur < pmax:
            iCount += 1
        else:
            break

    clump.minRadius = np.min(clump.radii)
    clump.maxRadius = np.max(clump.radii)
    clump.numSpheres = len(clump.radii)

    output = kwargs.get('output')
    if output is not None:
        np.savetxt(output, np.asarray(np.hstack((clump.positions, clump.radii))),
                   delimiter=",")  # In PyCharm this line seems to have an error but it does not. Known issue.

    outputVTK = kwargs.get('outputVTK')
    if outputVTK is not None:
        clump_to_VTK(clump, filepath=outputVTK)

    visualise = kwargs.get('visualise')
    if visualise is not None and visualise:
        clump_plotter_pyvista(clump)

    return mesh, clump
