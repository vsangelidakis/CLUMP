import numpy as np
import time
from scipy.spatial import Delaunay

"""
2017 Original author in MATLAB is Luigi Giaccari: https://www.mathworks.com/matlabcentral/fileexchange/63731-surface-reconstruction-from-scattered-points-cloud-open-surfaces?s_tid=srchtitle

2024 Translated from MATLAB to Python by A.U. Canbolat <utku.canbolat@fau.de>
"""

def MyRobustCrust(p):
    # error check
    if len(p.shape) > 2 or p.shape[1] != 3:
        raise ValueError("Input 3D points must be stored in a Nx3 array")

    # add points to the given ones, this is useful to create outside tetrahedrons
    start = time.time()
    p, nshield = AddShield(p)
#    print(f'Added Shield: {time.time() - start} s')

    start = time.time()
    # https://stackoverflow.com/questions/36604172/difference-between-matlab-delaunayn-and-scipy-delaunay
    tetr = matlab_delaunayn(p)

#    print(f'Delaunay Triangulation Time: {time.time() - start} s')

    # find triangles to tetrahedron and tetrahedron to triangles connectivity data
    start = time.time()
    t2tetr, tetr2t, t = Connectivity(tetr)
#    print(f'Connectivity Time: {time.time() - start} s')

    start = time.time()
    cc, r = CC(p, tetr)  # Circumcenters of tetrahedrons
#    print(f'Circumcenters Time: {time.time() - start} s')

    start = time.time()
    tbound, _ = Marking(p, tetr, tetr2t, t2tetr, cc, r, nshield)  # Flagging tetrahedrons as inside or outside
#    print(f'Walking Time: {time.time() - start} s')

    # reconstructed raw surface
    t = t[tbound]

    start = time.time()
    t, tnorm = ManifoldExtraction(t, p)
#    print(f'Manifold extraction Time: {time.time() - start} s')

    return t, tnorm


def TriAngle(p1, p2, p3, p4, planenorm):
    """
    Computes angle between two triangles.
    """

    # Determine if p4 is above or below the plane identified by planenorm and p3
    test = np.dot(planenorm, p4 - p3)

    # Compute triangle normals
    v21 = p1 - p2
    v31 = p3 - p1

    tnorm1 = np.array([
        v21[1] * v31[2] - v21[2] * v31[1],
        v21[2] * v31[0] - v21[0] * v31[2],
        v21[0] * v31[1] - v21[1] * v31[0]
    ])
    tnorm1 /= np.linalg.norm(tnorm1)

    v41 = p4 - p1
    tnorm2 = np.array([
        v21[1] * v41[2] - v21[2] * v41[1],
        v21[2] * v41[0] - v21[0] * v41[2],
        v21[0] * v41[1] - v21[1] * v41[0]
    ])
    tnorm2 /= np.linalg.norm(tnorm2)

    # Cosine of the angle
    alpha = np.dot(tnorm1, tnorm2)

    # Convert cosine value to angle
    alpha = np.arccos(np.clip(alpha, -1, 1))

    # Adjust the angle if p4 is below the plane
    if test < 0:
        alpha = alpha + 2 * (np.pi - alpha)

    # Check if we need to change the orientation of the second triangle
    testor = np.dot(planenorm, tnorm1)
    if testor > 0:
        tnorm2 = -tnorm2

    return alpha, tnorm2


def ManifoldExtraction(t, p):  # WORKS
    """
    :param t: A list of triangles. Each row in this list represents a triangle, with the three values being indices into the point list p.
    :param p: A list of 3D points.
    :return:    t: A modified list of triangles that represents the manifold surface.
                tnorm: Normals for the triangles in the manifold.
    """

    numt = t.shape[0]  # number of triangles
    vect = np.arange(numt)
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e, j = np.unique(np.sort(e, axis=1), axis=0, return_inverse=True)

    # Unique edges
    te = np.vstack([j[vect], j[vect + numt], j[vect + 2 * numt]]).T

    nume = e.shape[0]
    e2t = np.zeros((nume, 2), dtype=np.int32)
    ne = e.shape[0]
    # np_ = p.shape[0]
    count = np.zeros(ne, dtype=np.int32)
    etmapc = np.zeros((ne, 4), dtype=np.int32)

    for i in range(numt):
        i1, i2, i3 = te[i, :]
        etmapc[i1, count[i1]] = i
        etmapc[i2, count[i2]] = i
        etmapc[i3, count[i3]] = i

        count[i1] += 1
        count[i2] += 1
        count[i3] += 1

    etmap = [etmapc[i, :count[i]].tolist() for i in range(ne)]
    tkeep = np.zeros(numt, dtype=bool)

    tnorm = Tnorm(p, t)

    t1 = np.argmax(np.sum(p[t, 2], axis=1) / 3)
    if tnorm[t1, 2] < 0:
        tnorm[t1, :] = -tnorm[t1, :]

    tkeep[t1] = True
    efront = np.zeros(nume, dtype=np.int32)
    efront[:3] = te[t1, :]
    e2t[te[t1, :], 0] = t1
    nf = 2

    while nf > 0:
        k = efront[nf]
        if e2t[k, 1] > 0 or e2t[k, 0] < 1 or count[k] < 2:
            nf -= 1
            continue

        idtcandidate = etmap[k]
        t1 = e2t[k, 0]
        ttemp = t[t1, :]
        etemp = e[k, :]
        p1, p2 = etemp
        p3 = ttemp[(ttemp != p1) & (ttemp != p2)][0]

        alphamin = np.inf
        for i in range(len(idtcandidate)):
            t2 = idtcandidate[i]
            if t2 == t1:
                continue

            ttemp = t[t2, :]
            p4 = ttemp[(ttemp != p1) & (ttemp != p2)][0]
            alpha, tnorm2 = TriAngle(p[p1, :], p[p2, :], p[p3, :], p[p4, :], tnorm[t1, :])

            if alpha < alphamin:
                alphamin = alpha
                idt = t2
                tnorm[t2, :] = tnorm2

        tkeep[idt] = True
        for j in range(3):
            ide = te[idt, j]

            if ide == 0:
                break

            if e2t[ide, 0] < 1:  # Is it the first triangle for the current edge?
                efront[nf] = ide
                nf += 1
                e2t[ide, 0] = idt
            else:  # No, it is the second one
                efront[nf] = ide
                nf += 1
                e2t[ide, 1] = idt

        nf -= 1

    t = t[tkeep, :]
    tnorm = tnorm[tkeep, :]

    return t, tnorm


def CC(p, tetr):
    # Finds circumcenters from a set of tetrahedrons

    ntetr = tetr.shape[0]
    cc = np.zeros((ntetr, 3))

    # Points of tetrahedron
    p1 = p[tetr[:, 0]]
    p2 = p[tetr[:, 1]]
    p3 = p[tetr[:, 2]]
    p4 = p[tetr[:, 3]]

    # Vectors of tetrahedron edges
    v21 = p1 - p2
    v31 = p3 - p1
    v41 = p4 - p1

    # Solve the system using Cramer method
    d1 = np.sum(v41 * (p1 + p4) * 0.5, axis=1)
    d2 = np.sum(v21 * (p1 + p2) * 0.5, axis=1)
    d3 = np.sum(v31 * (p1 + p3) * 0.5, axis=1)

    det23 = v21[:, 1] * v31[:, 2] - v21[:, 2] * v31[:, 1]
    det13 = v21[:, 2] * v31[:, 0] - v21[:, 0] * v31[:, 2]
    det12 = v21[:, 0] * v31[:, 1] - v21[:, 1] * v31[:, 0]

    Det = v41[:, 0] * det23 + v41[:, 1] * det13 + v41[:, 2] * det12

    detx = d1 * det23 + v41[:, 1] * (-d2 * v31[:, 2] + v21[:, 2] * d3) + v41[:, 2] * (d2 * v31[:, 1] - v21[:, 1] * d3)

    dety = v41[:, 0] * (d2 * v31[:, 2] - v21[:, 2] * d3) + d1 * det13 + v41[:, 2] * (d3 * v21[:, 0] - v31[:, 0] * d2)

    detz = v41[:, 0] * (v21[:, 1] * d3 - d2 * v31[:, 1]) + v41[:, 1] * (d2 * v31[:, 0] - v21[:, 0] * d3) + d1 * det12

    cc[:, 0] = detx / Det
    cc[:, 1] = dety / Det
    cc[:, 2] = detz / Det

    # Circumradius
    r = np.sqrt(np.sum((p2 - cc) ** 2, axis=1))  # circumradius

    return cc, r


def Connectivity(tetr):
    # Gets connectivity relationships among tetrahedrons
    numt = tetr.shape[0]
    vect = np.arange(numt)
    t = np.vstack([tetr[:, [0, 1, 2]], tetr[:, [1, 2, 3]], tetr[:, [0, 2, 3]], tetr[:, [0, 1, 3]]])  # triangles not unique
    t, j = np.unique(np.sort(t, axis=1), return_inverse=True, axis=0)  # triangles
    t2tetr = np.vstack([j[vect], j[vect + numt], j[vect + 2 * numt], j[vect + 3 * numt]]).T  # each tetrahedron has 4 triangles

    # triang-to-tetr connectivity
    nume = t.shape[0]
    tetr2t = np.zeros((nume, 2), dtype=np.int32)
    count = np.ones(nume, dtype=np.int8)
    for k in range(numt):
        for j in range(4):
            ce = t2tetr[k, j]
            tetr2t[ce, count[ce] - 1] = k
            count[ce] += 1

    return t2tetr, tetr2t, t


def Marking(p, tetr, tetr2t, t2tetr, cc, r, nshield):
    # constants for the algorithm
    TOLLDIFF = .01
    INITTOLL = .99
    MAXLEVEL = 10 / TOLLDIFF
    BRUTELEVEL = MAXLEVEL - 50

    # preallocation
    np_ = p.shape[0] - nshield
    numtetr = tetr.shape[0]
    nt = tetr2t.shape[0]

    # First flag as outside tetrahedrons with Shield points
    deleted = np.any(tetr > np_, axis=1)
    checked = deleted.copy()
    onfront = np.zeros(nt, dtype=bool)

    for i in np.where(checked)[0]:
        onfront[t2tetr[i, :]] = True

    countchecked = np.sum(checked)

    # tolerances to mark as in or out
    toll = np.full(nt, INITTOLL)
    level = 0

    # intersection factor
    Ifact = IntersectionFactor(tetr2t, cc, r)

    ids = np.arange(nt)
    queue = ids[onfront]
    nt = len(queue)
    while countchecked < numtetr and level < MAXLEVEL:
        level += 1

        for i in range(nt):
            id_ = queue[i]
            tetr1, tetr2 = tetr2t[id_]

            if tetr2 == 0 or (checked[tetr1] and checked[tetr2]):
                onfront[id_] = False
                continue

            if Ifact[id_] >= toll[id_]:  # flag as equal
                if checked[tetr1]:
                    deleted[tetr2] = deleted[tetr1]
                    checked[tetr2] = True
                    countchecked += 1
                    onfront[t2tetr[tetr2]] = True
                else:
                    deleted[tetr1] = deleted[tetr2]
                    checked[tetr1] = True
                    countchecked += 1
                    onfront[t2tetr[tetr1]] = True
                onfront[id_] = False

            elif Ifact[id_] < -toll[id_]:  # flag as different
                if checked[tetr1]:
                    deleted[tetr2] = not deleted[tetr1]
                    checked[tetr2] = True
                    countchecked += 1
                    onfront[t2tetr[tetr2]] = True
                else:
                    deleted[tetr1] = not deleted[tetr2]
                    checked[tetr1] = True
                    countchecked += 1
                    onfront[t2tetr[tetr1]] = True
                onfront[id_] = False

            else:
                toll[id_] -= TOLLDIFF

        if level == BRUTELEVEL:
            print('Brute continuation necessary')
            onfront[np.any(t2tetr[~checked], axis=0)] = True

        queue = ids[onfront]
        nt = len(queue)

    # extract boundary triangles
    tbound = BoundTriangles(tetr2t, deleted)

    if level == MAXLEVEL:
        print(f'{level} th level was reached')
#    else:
#        print(f'{level} th level was reached')

#    print(f'{countchecked / numtetr * 100:.4f} % of tetrahedrons were checked')

    return tbound, Ifact


def AddShield(p):
    """
    Adds outside points to the given cloud forming outside tetrahedrons.

    Shield points are very good in detecting outside tetrahedrons. Unfortunately,
    Delaunay triangulation with these points can be even 50% slower.
    """

    # Find the bounding box
    maxx = np.max(p[:, 0])
    maxy = np.max(p[:, 1])
    maxz = np.max(p[:, 2])
    minx = np.min(p[:, 0])
    miny = np.min(p[:, 1])
    minz = np.min(p[:, 2])

    # Give offset to the bounding box
    step = max([maxx - minx, maxy - miny, maxz - minz])

    maxx += step
    maxy += step
    maxz += step
    minx -= step
    miny -= step
    minz -= step

    N = 10  # Number of points of the shield edge

    step /= (N * N)  # Decrease step, avoids not unique points

    nshield = N * N * 6

    # Creating a grid lying on the bounding box
    vx = np.linspace(minx, maxx, N)
    vy = np.linspace(miny, maxy, N)
    vz = np.linspace(minz, maxz, N)

    x, y = np.meshgrid(vx, vy)
    facez1 = np.column_stack([x.ravel(), y.ravel(), np.ones(N * N) * maxz])
    facez2 = np.column_stack([x.ravel(), y.ravel(), np.ones(N * N) * minz])

    x, y = np.meshgrid(vy, vz - step)
    facex1 = np.column_stack([np.ones(N * N) * maxx, x.ravel(), y.ravel()])
    facex2 = np.column_stack([np.ones(N * N) * minx, x.ravel(), y.ravel()])

    x, y = np.meshgrid(vx - step, vz)
    facey1 = np.column_stack([x.ravel(), np.ones(N * N) * maxy, y.ravel()])
    facey2 = np.column_stack([x.ravel(), np.ones(N * N) * miny, y.ravel()])

    # Add points to the p array
    pnew = np.vstack([p, facex1, facex2, facey1, facey2, facez1, facez2])

    return pnew, nshield


def BoundTriangles(tetr2t, deleted):
    """ Extracts boundary triangles from a set tetr2t connectivity and form the deleted vector which tells tetrahedrons that are marked as out
    """
    nt = tetr2t.shape[0]  # Number of total triangles

    tbound = np.ones((nt, 2), dtype=bool)  # Initialize to keep shape in next operation

    ind = tetr2t > 0  # Avoid null index
    tbound[ind] = deleted[tetr2t[ind]]  # Mark 1 for deleted 0 for kept tetrahedrons

    tbound = np.sum(tbound, axis=1) == 1  # Boundary triangles only have one tetrahedron

    return tbound


def IntersectionFactor(tetr2t, cc, r):
    nt = tetr2t.shape[0]
    Ifact = np.zeros((nt, 1))  # Intersection factor

    i = tetr2t[:, 1] > 0

    # Distance between circumcenters
    distcc = np.sum((cc[tetr2t[i, 0], :] - cc[tetr2t[i, 1], :]) ** 2, axis=1)

    # Intersection factor
    Ifact[i] = ((-distcc + r[tetr2t[i, 0]].flatten() ** 2 + r[tetr2t[i, 1]].flatten() ** 2) / (
            2 * r[tetr2t[i, 0]].flatten() * r[tetr2t[i, 1]].flatten()))[:, np.newaxis]

    return Ifact


def Tnorm(p, t):  # WORKS
    """ Computes normalized normals of triangles"""
    v21 = p[t[:, 0]] - p[t[:, 1]]
    v31 = p[t[:, 2]] - p[t[:, 0]]

    tnorm1 = np.zeros(t.shape)

    tnorm1[:, 0] = v21[:, 1] * v31[:, 2] - v21[:, 2] * v31[:, 1]
    tnorm1[:, 1] = v21[:, 2] * v31[:, 0] - v21[:, 0] * v31[:, 2]
    tnorm1[:, 2] = v21[:, 0] * v31[:, 1] - v21[:, 1] * v31[:, 0]

    L = np.sqrt(np.sum(tnorm1 ** 2, axis=1))

    tnorm1 = tnorm1 / L[:, None]

    return tnorm1


def merge_duplicate_points(X):
    """Merge out points that have coincident location"""
    dupes_found = False
    num_init_points = X.shape[0]

    # Get unique rows and their indices
    _, idx_map = np.unique(X, axis=0, return_index=True)

    num_unique_points = len(idx_map)
    if num_init_points > num_unique_points:
        # Undo the sort to preserve the ordering of points
        idx_map.sort()
        X = X[idx_map]
        dupes_found = True

    return X, dupes_found, idx_map


def matlab_delaunayn(x):
    if x is None:
        raise ValueError('Not Enough Inputs')

    n = x.shape[1]

    if n < 1:
        raise ValueError('X has Low Column Number')

    x, dupesfound, idxmap = merge_duplicate_points(x)

    if dupesfound:
        print('Warning: Duplicate Data Points')

    m, n = x.shape

    if m < n + 1:
        raise ValueError('Not Enough Points for Tessellation')

    if m == n + 1:
        t = np.arange(n + 1)

        # Enforce the orientation convention
        if n == 2 or n == 3:
            PredicateMat = np.ones((m, m))
            PredicateMat[:, :n] = x[t, :n]
            orient = np.linalg.det(PredicateMat)

            if n == 3:
                orient *= -1

            if orient < 0:
                t[1], t[2] = t[2], t[1]

        if dupesfound:
            t = idxmap[t]

        return t

    t = Delaunay(x, qhull_options="Qt Qbb Qc")  # Scipy's qhull
    t = t.simplices

    # Strip the zero volume simplices that may have been created by the presence of degeneracy
    mt, nt = t.shape
    v = np.ones(mt, dtype=bool)

    for i in range(mt):
        xa = x[t[i, :nt - 1]]
        xb = x[t[i, nt - 1]]

        val = np.abs(np.linalg.det(xa - xb))
        valtol = np.finfo(float).eps * np.max(np.abs(np.concatenate((xa.flatten(), xb.flatten()))))

        if val <= valtol:
            v[i] = False

    t = t[v]

    if dupesfound:
        t = idxmap[t]

    return t
