import numpy as np
import CLUMP.utils.RigidBodyParameters as RigidBodyParameters
import CLUMP.utils.STL_ReaderWriter as STL_ReaderWriter
from CLUMP.utils.ParticlePlotter import clump_plotter_pyvista
from CLUMP.utils.VTK_Writer import clump_to_VTK

"""
Implementation of the clump-generation concept proposed by Favier et al. (1999) [1]
2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

2021 MATLAB implementation by V. Angelidakis <v.angelidakis@qub.ac>
2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>

[1] Favier, J.F., Fard, M.H., Kremmer, M. and Raji, A.O., 1999. Engineering Computations: Int J for Computer-Aided Engineering, 16(4), pp.467-480.

The main concept of this methodology:

1. We import the geometry of a particle either as a surface mesh or a voxelated 3D image.
2. If a voxelated image is given, we transform it into a surface mesh, providing its vertices and faces (vertex connectivity).
3. We calculate the inertial characteristics of the particle and center it to its centroid and align it to its principal axes.
4. We create a number of N points along the longest particle axis and identify the particle vertices belonging to each of the newly formed (N+1) spans.
5. We generate a sphere centered to each of the N points. The radius of each sphere is calculated based on the distances of the vertices within the span of 
interest. 
The default behaviour considers the minimum distance, although an optional parameter (chooseDistance) exists that takes the values 'min' (default), 'avg' and 'max'.
"""


class Clump:
    def __init__(self):
        self.positions = np.empty((0, 3))
        self.radii = np.empty((0, 1))
        self.minRadius = None
        self.maxRadius = None
        self.numSpheres = None


def GenerateClump_Favier(inputGeom, N, **kwargs):
    """ Function to generate clumps using the concept proposed by Favier et al. (1999)
    :param inputGeom: Input geometry, given in one of the formats below:
                    1. Directory of .stl file (for surface meshes)
                    2. Directory of .mat file (for binary voxelated images)
                    3. Struct with fields {vertices,faces} (for surface meshes)
                    4. Struct with fields {img,voxel_size} (for binary voxelated images) where
                        - img:			[Nx x Ny x Nz] voxelated image
                        - voxel_size:	[1x3] voxel size in Cartesian space
    :param N: Number of spheres to be generated.
    :param kwargs: Can contain either of the optimal variables "chooseDistance", "output", "outputVTK", "visualise".
                - chooseDistance: Preferred method to specify the radii of the spheres, which can be either the minimum ('min'),
                the average ('avg') or the maximum ('max') distance of the vertices within the span of interest.
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

    # Configure input particle geometry based on the variable type of inputGeom
    # THIS PART NEEDS TO BE IMPLEMENTED. CURRENTLY ONLY ALLOW STL FILES

    clump = Clump()  # instentiate Clump object for later use

    F, P = STL_ReaderWriter.read_stl(inputGeom)  # read the STL file and get faces and vertices

    # Build "mesh" structure
    mesh = RigidBodyParameters.RBP(F, P)  # calculate the rigid body parameters based on F and P.

    ################################################################################################
    #                                   Main Body of the Function                                  #
    ################################################################################################

    # Center particle around its centroid and align it along its principal axes.
    rot = mesh.PAI
    # Singular value decomposition implementations in numpy and MATLAB give different results due to the
    # freedom of choosing the basis vectors. To make them same I added two lines below. They can be removed
    # without loss of generality. - Utku
    rot *= -1
    rot[:, 1] *= -1

    # Transform rotation matrix to align longest axis along X direction
    rot[:, [2, 0]] = rot[:, [0, 2]]

    P = P - mesh.centroid
    P = P @ rot

    X_extremas = np.array([np.min(P[:, 0]), np.max(P[:, 0])])
    a, b = X_extremas[0], X_extremas[1]

    nSegments = N

    endPoints = np.linspace(a, b, nSegments + 1)  # 4 endpoints for 3 segments
    start = endPoints[0:-1]  # 3 starting points
    stop = endPoints[1::]  # 3 stopping points
    midPoints = stop - (stop[0] - start[0]) / 2  # 3 middle points

    minDistance = np.zeros(midPoints.size)
    minDx = np.zeros(midPoints.size)

    p_struct = {}  # MATLAB's struct is equivalent to Python's dictionary
    for i in range(midPoints.size):  # For each midpoint
        temp_arr = np.array([])  # I initialize an emtpy array to stack the values required below

        # Find vertices within each sector
        for j in range(P.shape[0]):  # for each point of the particle surface (in principal axes)
            if endPoints[i] <= P[j, 0] <= endPoints[i + 1]:
                temp_arr = np.append(temp_arr, P[j, 0:3])  # I concatante the (3, 1) arrays
        p_struct[i] = temp_arr.reshape(temp_arr.size // 3, 3)  # then I reshape the temp_arr to add it to p_struct

        # Find closest distance of midpoint to any surface vertex
        xM, yM, zM = midPoints[i], 0, 0
        x1 = endPoints[0]
        x2 = endPoints[-1]

        # Closest distance between midpoint and particle surface
        minDistance[i] = np.min(
            np.sqrt(
                np.square(P[:, 0] - xM) + np.square(P[:, 1] - yM) + np.square(P[:, 2] - zM)))

        # Closest distane between midpoint and particle X limits
        minDx[i] = np.min((np.abs(x1 - xM), np.abs(x2 - xM)))

    # Build "clump" structure
    radius = np.zeros(midPoints.size)
    for i in range(midPoints.size):  # I skipped the isempty query.
        distance = np.sqrt(np.square(p_struct[i][:, 0] - midPoints[i]) + np.square(p_struct[i][:, 1]) + np.square(
            p_struct[i][:, 2]))

        # I wrote the switch-case expression in old-school fashion for backward compability
        chooseDistance = kwargs.get('chooseDistance')
        if chooseDistance == "min":
            radius[i] = np.min(distance)
        elif chooseDistance is None:
            radius[i] = np.min(distance)
        elif chooseDistance == "avg":
            radius[i] = np.mean(distance)
        elif chooseDistance == "max":
            radius[i] = np.max(distance)
        else:
            print("Wrong optional parameter type for chooseDistance.")

        radius[i] = np.min((radius[i], minDx[i]))

        clump.positions = np.vstack((clump.positions, np.array([midPoints[i], 0, 0])))
        clump.radii = np.vstack((clump.radii, radius[i]))

    # Transform the mesh and the clump coordinates back to the initial (non-principal) system
    P = P @ np.transpose(rot)
    P += mesh.centroid

    clump.positions = clump.positions @ np.transpose(rot)
    clump.positions += mesh.centroid

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
