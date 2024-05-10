import numpy as np

"""
This script calculates the rigid body parameters of objects imported as surface meshes. It provides their volume, centroid, inertia tensor, principal inertia tensor, principal inertia vectors, moments.

2021 Original author in MATLAB is Anton Semechko: https://www.mathworks.com/matlabcentral/fileexchange/48913-rigid-body-parameters-of-closed-surface-meshes

2024 Translated from MATLAB to Python by U.A. Canbolat <utku.canbolat@fau.de>
"""

# this is an ugly design. I will fix it later
class RBP:
    def __init__(self, F, P):
        self.F, self.P = F, P

        results = self.calculate_RBP()

        self.moments = results[0:10]
        self.volume = results[0]
        self.centroid = np.array(results[1:4])/results[0]
        self.eigs = results[-1]
        self.PAI = results[-2]
        self.inertia_tensor = results[-3]

    def calculate_RBP(self):
        """ Area weighted face normals (magnitude = 2 x triangle area)"""
        X1 = self.P[self.F[:, 0], :]
        X2 = self.P[self.F[:, 1], :]
        X3 = self.P[self.F[:, 2], :]

        FN = np.cross(X2 - X1, X3 - X1)

        # Zeroth order moment (same as volume enclosed by the mesh) ---------------
        C = (X1 + X2 + X3) / 3.0
        m000 = np.sum(np.sum(np.multiply(FN, C))) / 6.0

        epsilon = np.finfo(float).eps  # Machine epsilon
        if m000 < epsilon:
            print(
                'Mesh has negative volume, indicating that face normal orientations are reversed (i.e., pointing into the surface).')

        # % First order moments (together they specify the centroid of the region
        # % enclosed by the mesh) ---------------------------------------------------

        x1, y1, z1 = X1[:, 0], X1[:, 1], X1[:, 2]
        x2, y2, z2 = X2[:, 0], X2[:, 1], X2[:, 2]
        x3, y3, z3 = X3[:, 0], X3[:, 1], X3[:, 2]

        x_2 = (np.multiply((x1 + x2), (x2 + x3)) + np.square(x1) + np.square(x3)) / 12.0
        y_2 = (np.multiply((y1 + y2), (y2 + y3)) + np.square(y1) + np.square(y3)) / 12.0
        z_2 = (np.multiply((z1 + z2), (z2 + z3)) + np.square(z1) + np.square(z3)) / 12.0

        xy = (np.multiply((x1 + x2 + x3), (y1 + y2 + y3)) + np.multiply(x1, y1) + np.multiply(x2, y2) + np.multiply(x3,
                                                                                                                    y3)) / 24.0
        xz = (np.multiply((x1 + x2 + x3), (z1 + z2 + z3)) + np.multiply(x1, z1) + np.multiply(x2, z2) + np.multiply(x3,
                                                                                                                    z3)) / 24.0
        yz = (np.multiply((y1 + y2 + y3), (z1 + z2 + z3)) + np.multiply(y1, z1) + np.multiply(y2, z2) + np.multiply(y3,
                                                                                                                    z3)) / 24.0

        m100 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((x_2, 2 * xy, 2 * xz)))))) / 6.0
        m010 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((2 * xy, y_2, 2 * yz)))))) / 6.0
        m001 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((2 * xz, 2 * yz, z_2)))))) / 6.0

        # % Second order moments (used to determine elements of the inertia tensor)
        # % -------------------------------------------------------------------------

        x_3 = (np.multiply((x1 + x2 + x3), (np.square(x1) + np.square(x2) + np.square(x3))) + np.multiply(x1,
                                                                                                          np.multiply(
                                                                                                              x2,
                                                                                                              x3))) / 20.0
        y_3 = (np.multiply((y1 + y2 + y3), (np.square(y1) + np.square(y2) + np.square(y3))) + np.multiply(y1,
                                                                                                          np.multiply(
                                                                                                              y2,
                                                                                                              y3))) / 20.0
        z_3 = (np.multiply((z1 + z2 + z3), (np.square(z1) + np.square(z2) + np.square(z3))) + np.multiply(z1,
                                                                                                          np.multiply(
                                                                                                              z2,
                                                                                                              z3))) / 20.0

        x_2y = (np.multiply((3 * y1 + y2 + y3), np.square(x1)) + np.multiply((y1 + 3 * y2 + y3),
                                                                             np.square(x2)) + np.multiply(
            (y1 + y2 + 3 * y3), np.square(x3))
                + np.multiply((2 * y1 + 2 * y2 + y3), np.multiply(x1, x2))
                + np.multiply((2 * y1 + y2 + 2 * y3), np.multiply(x1, x3))
                + np.multiply((y1 + 2 * y2 + 2 * y3), np.multiply(x2, x3))) / 60.0

        x_2z = (np.multiply((3 * z1 + z2 + z3), np.square(x1)) + np.multiply((z1 + 3 * z2 + z3),
                                                                             np.square(x2)) + np.multiply(
            (z1 + z2 + 3 * z3), np.square(x3))
                + np.multiply((2 * z1 + 2 * z2 + z3), np.multiply(x1, x2))
                + np.multiply((2 * z1 + z2 + 2 * z3), np.multiply(x1, x3))
                + np.multiply((z1 + 2 * z2 + 2 * z3), np.multiply(x2, x3))) / 60.0

        y_2x = (np.multiply((3 * x1 + x2 + x3), np.square(y1)) + np.multiply((x1 + 3 * x2 + x3),
                                                                             np.square(y2)) + np.multiply(
            (x1 + x2 + 3 * x3), np.square(y3))
                + np.multiply((2 * x1 + 2 * x2 + x3), np.multiply(y1, y2))
                + np.multiply((2 * x1 + x2 + 2 * x3), np.multiply(y1, y3))
                + np.multiply((x1 + 2 * x2 + 2 * x3), np.multiply(y2, y3))) / 60.0

        y_2z = (np.multiply((3 * z1 + z2 + z3), np.square(y1)) + np.multiply((z1 + 3 * z2 + z3),
                                                                             np.square(y2)) + np.multiply(
            (z1 + z2 + 3 * z3), np.square(y3))
                + np.multiply((2 * z1 + 2 * z2 + z3), np.multiply(y1, y2))
                + np.multiply((2 * z1 + z2 + 2 * z3), np.multiply(y1, y3))
                + np.multiply((z1 + 2 * z2 + 2 * z3), np.multiply(y2, y3))) / 60.0

        z_2y = (np.multiply((3 * y1 + y2 + y3), np.square(z1)) + np.multiply((y1 + 3 * y2 + y3),
                                                                             np.square(z2)) + np.multiply(
            (y1 + y2 + 3 * y3), np.square(z3))
                + np.multiply((2 * y1 + 2 * y2 + y3), np.multiply(z1, z2))
                + np.multiply((2 * y1 + y2 + 2 * y3), np.multiply(z1, z3))
                + np.multiply((y1 + 2 * y2 + 2 * y3), np.multiply(z2, z3))) / 60.0

        z_2x = (np.multiply((3 * x1 + x2 + x3), np.square(z1)) + np.multiply((x1 + 3 * x2 + x3),
                                                                             np.square(z2)) + np.multiply(
            (x1 + x2 + 3 * x3), np.square(z3))
                + np.multiply((2 * x1 + 2 * x2 + x3), np.multiply(z1, z2))
                + np.multiply((2 * x1 + x2 + 2 * x3), np.multiply(z1, z3))
                + np.multiply((x1 + 2 * x2 + 2 * x3), np.multiply(z2, z3))) / 60.0

        xyz = (np.multiply((x1 + x2 + x3), np.multiply((y1 + y2 + y3), (z1 + z2 + z3)))
               - 0.5 * np.multiply((np.multiply(y2, z3) + np.multiply(y3, z2) - 4 * np.multiply(y1, z1)), x1)
               - 0.5 * np.multiply((np.multiply(y1, z3) + np.multiply(y3, z1) - 4 * np.multiply(y2, z2)), x2)
               - 0.5 * np.multiply((np.multiply(y1, z2) + np.multiply(y2, z1) - 4 * np.multiply(y3, z3)), x3)) / 60.0

        m110 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((x_2y, y_2x, 2 * xyz)))))) / 6.0
        m101 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((x_2z, 2 * xyz, z_2x)))))) / 6.0
        m011 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((2 * xyz, y_2z, z_2y)))))) / 6.0

        m200 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((x_3, 3 * x_2y, 3 * x_2z)))))) / 9.0
        m020 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((3 * y_2x, y_3, 3 * y_2z)))))) / 9.0
        m002 = np.sum(np.sum(np.multiply(FN, np.transpose(np.vstack((3 * z_2x, 3 * z_2y, z_3)))))) / 9.0

        # % Inertia tensor ----------------------------------------------------------

        Ixx = m020 + m002 - (m010 ** 2 + m001 ** 2) / m000
        Iyy = m200 + m002 - (m100 ** 2 + m001 ** 2) / m000
        Izz = m200 + m020 - (m100 ** 2 + m010 ** 2) / m000
        Ixy = m110 - m100 * m010 / m000
        Ixz = m101 - m100 * m001 / m000
        Iyz = m011 - m010 * m001 / m000

        Inertia = np.array([[Ixx, -Ixy, -Ixz],
                            [-Ixy, Iyy, -Iyz],
                            [-Ixz, -Iyz, Izz]])

        # ------------------------------------------------------------------------------

        V, D, _ = np.linalg.svd(Inertia)

        return m000, m100, m010, m001, m110, m101, m011, m200, m020, m002, Inertia, V, D


