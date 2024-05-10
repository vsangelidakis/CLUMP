import numpy as np
from stl.mesh import Mesh

"2024 Python implementation by A.U. Canbolat <utku.canbolat@fau.de>"

def read_stl(stl_dir, is_duplicated=False):
    """ Read stl files into mesh. mesh contains all the required information"""
    mesh = Mesh.from_file(stl_dir)

    # extract the vectors of vertices
    P = mesh.points.reshape(-1, 3)

    # face enumeration
    F = np.arange(0, P.shape[0], dtype=np.int64).reshape(mesh.points.shape[0], 3)

    if is_duplicated:
        return F, P
    else:
        # now take unique values to do stl.SlimVerts.m's job
        # it discards the common vertices to avoid duplication.
        P_unique, indices = np.unique(P, return_inverse=True, axis=0)

        F_unique = indices.reshape(-1, 3)

        return F_unique, P_unique

def compute_normal(v1, v2, v3):
    """ Compute the normal vector of a triangle given its vertices.

    Parameters:
        v1, v2, v3 (tuple): Vertices of the triangle.

    Returns:
        tuple: Normal vector (nx, ny, nz).
    """
    a = np.array(v2) - np.array(v1)
    b = np.array(v3) - np.array(v1)
    normal = np.cross(a, b)
    return normal / np.linalg.norm(normal)

def write_stl(filename, vertices, faces):
    """ Write vertices and faces to an STL file using stl.mesh.Mesh.

    Parameters:
        filename (str): Name of the STL file.
        vertices (list): List of 3D vertices [(x1, y1, z1), (x2, y2, z2), ...].
        faces (list): List of faces [(v1, v2, v3), ...], where v1, v2, v3 are indices of vertices.
    """
    # Create the mesh
    triangles = []
    for face in faces:
        v1, v2, v3 = face
        normal = compute_normal(vertices[v1], vertices[v2], vertices[v3])
        triangle = np.array([vertices[v1], vertices[v2], vertices[v3]])
        triangles.append((normal, triangle))
    
    # Create the mesh object
    stl_mesh = Mesh(np.zeros(len(triangles), dtype=Mesh.dtype))
    for i, (normal, triangle) in enumerate(triangles):
        stl_mesh.vectors[i] = triangle
        stl_mesh.normals[i] = normal
    
    # Write to file
    stl_mesh.save(filename)

# Example usage:
if __name__ == "__main__":
    # Read STL file to get vertices and faces
    F, P = read_stl("input.stl")

    # Write to new STL file
    write_stl("output.stl", P.tolist(), F.tolist())


