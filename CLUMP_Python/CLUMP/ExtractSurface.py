import numpy as np
from scipy.spatial import ConvexHull
from CLUMP.utils.MyRobustCrust import MyRobustCrust
from CLUMP.utils.ParticlePlotter import mesh_plotter_trimesh

"""
Tesselation of the surface of a clump into a surface mesh
2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

2021 MATLAB implementation by V. Angelidakis <v.angelidakis@qub.ac>
2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>
"""


def sphereContact(sphere1, sphere2):
    """ Function to perform contact detection between two spheres
    inContact: boolean: whether sphere1 and sphere2 intersect
    :param sphere1: [1 x 4] [x,y,z,r]:	test sphere 1
    :param sphere2: [1 x 4] [x,y,z,r]:	test sphere 2
    :return: boolean
    """
    d0 = np.linalg.norm(sphere2[0:3] - sphere1[0:3])  # Centroidal distance of the spheres
    if sphere1[3] + sphere2[3] > d0 > np.abs(sphere1[3] - sphere2[3]):
        return True
    else:
        return False


def makeSphere(X, Y, Z, radius, N):
    """ Function to create a surface mesh of a sphere with radius r, centered at (x,y,z) with N vertices
    :param X:
    :param Y:
    :param Z:
    :param radius:
    :param N:
    :return: vertices/faces
    """
    vertices = np.zeros((N, 3))
    inc = np.pi * (3 - np.sqrt(5))
    off = 2 / N
    for k in range(N):
        y = k * off - 1 + off / 2
        r = np.sqrt(1 - y ** 2)
        phi = k * inc
        vertices[k, 0:3] = np.array([np.cos(phi) * r * radius, y * radius, np.sin(phi) * r * radius])

    vertices += [X, Y, Z]

    # Compute the convex hull using scipy
    hull = ConvexHull(vertices, qhull_options="Qt")  # https://stackoverflow.com/questions/33615249/convexhull-orders-differently-in-matlab-and-python

    # 1. Sort the faces based on vertex indices
    sorted_faces = np.sort(hull.simplices, axis=1)
    sorted_indices = np.lexsort((sorted_faces[:, 2], sorted_faces[:, 1], sorted_faces[:, 0]))
    sorted_faces = sorted_faces[sorted_indices]

    # 2. Ensure consistent orientation
    for i in range(len(sorted_faces)):
        normal = np.dot(hull.equations[i, :3], hull.points[sorted_faces[i]].mean(axis=0))
        if normal > 0:  # Inconsistent orientation
            sorted_faces[i] = sorted_faces[i][[0, 2, 1]]  # Swap the last two indices

    return vertices, sorted_faces


def spherePotential(point, sphere, allowZero):
    """ Function to determine whether a point is inside a sphere
    :param point: [1 x 3] x,y,z:	test point
    :param sphere: [1 x 4] x,y,z,r	sphere of interest
    :param allowZero: boolean: whether to consider 0 values as contact, i.e. returning true
    :return: isInside: boolean: whether the test point is inside the sphere of interest
    """
    if allowZero:
        isInside = np.sqrt(((sphere[0] - point[0]) ** 2 + (sphere[1] - point[1]) ** 2 + (sphere[2] - point[2]) ** 2) / (sphere[3] ** 2)) - 1 <= 0
    else:
        isInside = np.sqrt(((sphere[0] - point[0]) ** 2 + (sphere[1] - point[1]) ** 2 + (sphere[2] - point[2]) ** 2) / (sphere[3] ** 2)) - 1 < 0

    return isInside


def ExtractSurface(clump, N_sphere, N_circle, visualise):
    """ Function to tessellate the surface of a clump into a surface mesh
    :param clump: either "clump" object or N x 4 matrix with columns of [x,y,z,r], where x,y,z the centroid of each sphere and r its radius
    :param N_sphere: Number of vertices on the surface of each member-sphere of the clump
    :param N_circle: Number of vertices on the circle defined as the intersection of two overlapping spheres
    :param visualise: Boolean whether to plot the generated surface mesh of the clump surface
    :return: faces: faces of generated surface mesh
             vertices : vertices of generated surface mesh
    """

    # Check format of input
    if isinstance(clump, dict):
        if 'positions' in clump and 'radii' in clump:
            if clump['positions'].shape[1] < 3:
                raise ValueError('Invalid format; clump.positions should have size N x 3.')
            if clump['radii'].shape[1] > 1:
                raise ValueError('Invalid format; clump.radii should have size N x 1.')
            spheresList = np.hstack((clump['positions'], clump['radii']))
        else:
            raise ValueError('Invalid format! The dictionary should have keys "positions" and "radii".')
    elif isinstance(clump, np.ndarray):
        if clump.shape[1] != 4:
            raise ValueError('Invalid format, should be x,y,z,r.')
        spheresList = clump
    else:
        raise ValueError('Invalid format for clump.')

    x, y, z, r = spheresList[:, 0], spheresList[:, 1], spheresList[:, 2], spheresList[:, 3]

    # Contact detection between all spheres (all possible combinations) - Record interactions
    interactions = []
    for i in range(spheresList.shape[0] - 1):
        for j in range(i + 1, spheresList.shape[0]):
            if i == j:
                continue
            inContact = sphereContact(spheresList[i], spheresList[j])
            if inContact:
                interactions.append([i, j])

    # Generate points for each sphere
    S = []
    for i in range(spheresList.shape[0]):
        vertices, faces = makeSphere(x[i], y[i], z[i], r[i], N_sphere)
        S.append({'vertices': vertices, 'faces': faces})

    # Calculate points on the intersection of each pair of interacting spheres
    for i in range(len(interactions)):
        n = spheresList[interactions[i][1], :3] - spheresList[interactions[i][0], :3]
        d = np.linalg.norm(n)
        n /= d

        r1 = spheresList[interactions[i][0], 3]
        r2 = spheresList[interactions[i][1], 3]

        h = np.sqrt((2 * r1 * d)**2 - (r1**2 + d**2 - r2**2)**2) / (2 * d)
        alph = np.arccos((r1**2 + d**2 - r2**2) / (2 * r1 * d))
        h1 = r1 * (1 - np.cos(alph))

        C = spheresList[interactions[i][0], :3] + n * (r1 - h1)

        n3 = n
        n1 = np.array([n[2], 0, -n[0]])  # Vector perpendicular to n
        if np.linalg.norm(n1) == 0:
            n1 = np.array([n[1], 0, -n[0]])
        n1 /= np.linalg.norm(n1)
        n2 = np.cross(n3, n1)

        a = np.linspace(-np.pi, np.pi, N_circle // 2)
        px = C[0] + h * (n1[0] * np.cos(a) + n2[0] * np.sin(a))
        py = C[1] + h * (n1[1] * np.cos(a) + n2[1] * np.sin(a))
        pz = C[2] + h * (n1[2] * np.cos(a) + n2[2] * np.sin(a))

        circlevertices = np.vstack((px, py, pz)).T

        S[interactions[i][0]]['circlevertices'] = circlevertices
        S[interactions[i][0]]['vertices'] = np.vstack((S[interactions[i][0]]['vertices'], circlevertices))

    # Perform contact detection to detect and delete points of each sphere
    # that are included in other spheres, in order to get only the points
    # on the surface of the clump
    for i in range(len(interactions)):
        sphere1 = S[interactions[i][0]]
        sphere2 = S[interactions[i][1]]

        # For interaction [sphere1,sphere2], check which vertices of sphere1 are inside sphere2
        inside_indices = []
        for j in range(sphere1['vertices'].shape[0]):
            if spherePotential(sphere1['vertices'][j], spheresList[interactions[i][1]], True):
                inside_indices.append(j)
        sphere1['vertices'] = np.delete(sphere1['vertices'], inside_indices, axis=0)

        # For interaction [sphere1,sphere2], check which vertices of sphere2 are inside sphere1
        inside_indices = []
        for j in range(sphere2['vertices'].shape[0]):
            if spherePotential(sphere2['vertices'][j], spheresList[interactions[i][0]], True):
                inside_indices.append(j)
        sphere2['vertices'] = np.delete(sphere2['vertices'], inside_indices, axis=0)

    # Collect vertices from all spheres in one variable
    vertices = np.vstack([s['vertices'] for s in S])
    vertices = np.unique(vertices, axis=0)

    # Generate mesh using the Crust algorithm (Amenta et al, 1999)
    faces, _ = MyRobustCrust(vertices)

    if visualise:
         mesh_plotter_trimesh(vertices, faces, spheresList)
    
    return faces, vertices
