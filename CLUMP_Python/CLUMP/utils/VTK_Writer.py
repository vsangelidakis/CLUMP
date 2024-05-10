import numpy as np
import os

"2024 Python implementation by A.U. Canbolat <utku.canbolat@fau.de>"

def clump_to_VTK(clump, filepath=os.path.join(os.getcwd(), "spheres.vtk")):
    """ Write sphere data to a VTK file without using external libraries.

    Parameters:
    spheres (np.array): Numpy array containing sphere data (positions and radii).
    output_dir (str): Directory where the VTK file will be saved.
    filename (str): Name of the VTK file.
    """

    # VTK header
    header = "# vtk DataFile Version 3.0\n"
    title = "Sphere data\n"
    datatype = "ASCII\n"
    structure = "DATASET POLYDATA\n"

    clump = np.hstack((clump.positions, clump.radii))

    # Count total number of spheres
    number_of_spheres = clump.shape[0]

    # Writing points
    points = f"POINTS {number_of_spheres} float\n"
    for sphere in clump:
        points += f"{sphere[0]} {sphere[1]} {sphere[2]}\n"

    # Writing point data (radii as scalars)
    point_data = f"POINT_DATA {number_of_spheres}\n"
    scalars = "SCALARS radius float 1\nLOOKUP_TABLE default\n"
    for sphere in clump:
        scalars += f"{sphere[3]}\n"

    # Combine sections
    vtk_content = header + title + datatype + structure + points + point_data + scalars

    # Write to file
    with open(filepath, 'w') as file:
        file.write(vtk_content)
