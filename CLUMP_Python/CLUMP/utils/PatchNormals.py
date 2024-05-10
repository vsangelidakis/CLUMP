import numpy as np
from numpy.matlib import repmat

"""
2017 Original author in MATLAB is Dirk-Jan Kroon (2024). Patch Normals (https://www.mathworks.com/matlabcentral/fileexchange/24330-patch-normals), MATLAB Central File Exchange. Retrieved April 25, 2024. 

2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>
"""

def patch_normals(F, P):
    """ Function to calculate the normal vectors of faces for a surface mesh."""
    Fa = F[:, 0]
    Fb = F[:, 1]
    Fc = F[:, 2]

    e1 = P[Fa, :] - P[Fb, :]
    e2 = P[Fb, :] - P[Fc, :]
    e3 = P[Fc, :] - P[Fa, :]

    e1_norm = np.divide(e1, np.transpose(
        repmat(np.sqrt(np.square(e1[:, 0]) + np.square(e1[:, 1]) + np.square(e1[:, 2])), 3, 1)))
    e2_norm = np.divide(e2, np.transpose(
        repmat(np.sqrt(np.square(e2[:, 0]) + np.square(e2[:, 1]) + np.square(e2[:, 2])), 3, 1)))
    e3_norm = np.divide(e3, np.transpose(
        repmat(np.sqrt(np.square(e3[:, 0]) + np.square(e3[:, 1]) + np.square(e3[:, 2])), 3, 1)))

    def elementwise_dot(mat1, mat2):  # I MUST FIND A WAY TO DO THIS OPERATION FASTER. NEED TO BE VECTORIZED.
        # mat1 and mat2 must be in same shape
        temp_mat = np.zeros(mat1.shape[0])
        for j in range(mat1.shape[0]):
            temp_mat[j] = np.dot(mat1[j, :], mat2[j, :])
        return temp_mat

    angle = np.transpose(np.array([np.arccos(elementwise_dot(e1_norm, -e3_norm)),
                                   np.arccos(elementwise_dot(e2_norm, -e1_norm)),
                                   np.arccos(elementwise_dot(e3_norm, -e2_norm))]))
    normal = np.cross(e1, e3)

    vertice_normals = np.zeros((P.shape[0], 3))

    for i in range(Fa.shape[0]):
        vertice_normals[Fa[i], :] += normal[i, :] * angle[i, 0]
        vertice_normals[Fb[i], :] += normal[i, :] * angle[i, 1]
        vertice_normals[Fc[i], :] += normal[i, :] * angle[i, 2]

    epsilon = np.finfo(float).eps  # Machine epsilon
    V_norm = np.sqrt(
        np.square(vertice_normals[:, 0]) + np.square(vertice_normals[:, 1]) + np.square(
            vertice_normals[:, 2])) + epsilon

    Nx = vertice_normals[:, 0] / V_norm
    Ny = vertice_normals[:, 1] / V_norm
    Nz = vertice_normals[:, 2] / V_norm

    N = np.transpose(np.vstack((Nx, np.vstack((Ny, Nz)))))

    return N
