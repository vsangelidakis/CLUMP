import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import trimesh

"2024 Python implementation by A.U. Canbolat and V. Angelidakis  <utku.canbolat@fau.de>"


def clump_plotter(P, F, clump):
    """Plots 3D spheres using trisurf"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create patch
    ax.plot_trisurf(P[:, 0], P[:, 1], P[:, 2], triangles=F, color=[0, 1, 0, 0.2], edgecolor='none')

    # Set properties
    ax.set_box_aspect([np.ptp(a) for a in [P[:, 0], P[:, 1], P[:, 2]]])  # axis equal
    ax.grid(True)  # grid on

    # Plot spheres
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]  # sphere mesh
    for i in range(len(clump.radii)):
        x, y, z, r = clump.positions[i, 0], clump.positions[i, 1], clump.positions[i, 2], clump.radii[i]
        X = r * np.cos(u) * np.sin(v) + x
        Y = r * np.sin(u) * np.sin(v) + y
        Z = r * np.cos(v) + z
        ax.plot_surface(X, Y, Z, color=np.random.rand(3), edgecolor='none')

    plt.show()


def add_sphere_to_plotter(plotter, center, radius, color, opacity, phi_res, theta_res):
    """Utility function to add a sphere to the plotter."""
    sphere = pv.Sphere(radius=radius, center=center, phi_resolution=phi_res, theta_resolution=theta_res)
    plotter.add_mesh(sphere, color=color, opacity=opacity)


def clump_plotter_pyvista(clump, opacity=1.0, phi_res=50, theta_res=50):
    """Plots 3D spheres using PyVista with increased opacity, a white background, and a finer mesh.
    Each sphere is defined by its center (x, y, z) and radius r.
    """
    spheres_data = np.hstack((clump.positions, clump.radii))
    plotter = pv.Plotter()

    colors = [plt.cm.jet(i / len(spheres_data)) for i in range(len(spheres_data))]  # Color mapping for spheres

    for sphere_data, color in zip(spheres_data, colors):
        add_sphere_to_plotter(plotter, sphere_data[:3], sphere_data[3], color, opacity, phi_res, theta_res)

    plotter.show(interactive=True)


def mesh_plotter_trimesh(vertices, faces, spheresList, phi_res=50, theta_res=50):
    """Plots a wireframe mesh and spheres using trimesh and PyVista."""
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    pv_mesh = pv.wrap(mesh)

    plotter = pv.Plotter()
    plotter.add_mesh(pv_mesh, style='wireframe', color='black', line_width=2)

    for (x, y, z, r) in spheresList:
        add_sphere_to_plotter(plotter, (x, y, z), r, 'white', 0.5, phi_res, theta_res)

    plotter.camera_position = [
        (np.min(vertices[:, 0]) + np.max(vertices[:, 0])) / 2,
        (np.min(vertices[:, 1]) + np.max(vertices[:, 1])) / 2,
        np.max(vertices[:, 2]) * 2
    ]

    plotter.show()
