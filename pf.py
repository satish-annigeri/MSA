import numpy as np
import numpy.typing as npt
from scipy.linalg import solve

Matrix = npt.NDArray[np.float_]


# Utility functions. Not directly related to DSM
def print_mat(header: str, k: Matrix) -> None:
    """
    Print the array k, preceded by a header string. Each element of k is
    printed using the format string fmt
    """
    m = k.shape[0]
    n = k.shape[1]
    print(f"{header} Size: {m} x {n}")
    for i in range(m):
        for j in range(n):
            print(f"{k[i, j]:12.4f}", end="")
        print()


def atba(a, b):
    return np.dot(np.dot(a.T, b), a)


# Location Matrix. Related to structure
def pf_calclm(n, bc):
    """
    Calculate location matrix containing degree of freedom numbers of each node
    """
    lm = np.zeros((n, 3), dtype=int)
    nd = 0
    ns = bc.shape[0]
    for i in range(ns):
        node = bc[i, 0]
        lm[node - 1, 0:3] = bc[i, 1:4]

    for node in range(n):
        for j in range(3):
            if lm[node, j] == 0:
                nd += 1
                lm[node, j] = nd
            else:
                lm[node, j] = 0
    return lm, nd


# Stiffness matrix in local coordinate system. Related to element
def pf_stiff(E: float, A: float, Iz: float, L: float) -> Matrix:
    """
    Calculate stiffness of a plane frame member in local coordinate system
    """
    k = np.zeros((6, 6), dtype=float)
    k[0, 0] = float(E) * A / L
    k[0, 3] = -k[0, 0]

    k[1, 1] = (12.0 * E * Iz) / L**3
    k[1, 2] = (6.0 * E * Iz) / L**2
    k[1, 4] = -k[1, 1]
    k[1, 5] = k[1, 2]

    k[2, 2] = (4.0 * E * Iz) / L
    k[2, 4] = -k[1, 2]
    k[2, 5] = k[2, 2] / 2.0

    k[3, 3] = k[0, 0]

    k[4, 4] = k[1, 1]
    k[4, 5] = -k[1, 2]

    k[5, 5] = k[2, 2]

    for i in range(6):
        for j in range(i):
            k[i, j] = k[j, i]

    return k


# Rotation matrix. Related to element
def pf_calcrot(dc):
    """
    Calculate rotation matrix of a plane frame member
    """
    r = np.zeros((6, 6), dtype=float)
    cx = dc[0]
    cy = dc[1]
    r[0, 0] = cx
    r[0, 1] = cy
    r[1, 0] = -cy
    r[1, 1] = cx
    r[2, 2] = 1
    r[3:6, 3:6] = r[0:3, 0:3]
    return r


# Length and direction cosines of a member. Related to element
def pf_calclen(xy1, xy2):
    """
    Calculate length and direction cosines of a plane frame member
    """
    delta = xy2 - xy1
    L = np.sqrt(np.sum(delta**2))
    dc = delta / L
    return L, dc


# Utility function to get number of end joints. Related to element
def pf_get_endjts(imem, conn):
    """Return numbers of first and second end of a plane frame member"""
    jt1 = conn[imem - 1, 0]
    jt2 = conn[imem - 1, 1]
    return jt1, jt2


# Utility function to get coordinates of end joints. Related to element
def pf_get_endcoord(imem, xy, conn):
    """Return the coordinates of the first and second end of a plane frame member"""
    jt1, jt2 = pf_get_endjts(imem, conn)
    xy1 = xy[jt1 - 1, :]
    xy2 = xy[jt2 - 1, :]
    return xy1, xy2


# Utility function to get member properties. Related to element
def pf_get_memprop(imem, xy, conn, mprop):
    """Return the properties of specified plane frame member"""
    m = conn[imem - 1, 2] - 1
    E = mprop[m, 0]
    A = mprop[m, 1]
    Iz = mprop[m, 2]
    xy1, xy2 = pf_get_endcoord(imem, xy, conn)
    L, dc = pf_calclen(xy1, xy2)
    return E, A, Iz, L, dc


# Utility function to get dof numbers at end joints. Related to element
def pf_get_dof(imem, conn, lm):
    """Return degree of freedom numbers of the ends of a plane frame member"""
    jt1, jt2 = pf_get_endjts(imem, conn)
    memdof = np.array([0, 0, 0, 0, 0, 0])
    memdof[0:3] = lm[jt1 - 1, :]
    memdof[3:6] = lm[jt2 - 1, :]
    return memdof


# Stiffness matrix in structure coordinate system. Related to element
def pf_gstiff(imem, xy, conn, mprop):
    """
    Calculate the stiffness matrix of a plane frame member in strucutre
    coordinate system
    """
    E, A, Iz, L, dc = pf_get_memprop(imem, xy, conn, mprop)
    r = pf_calcrot(dc)
    k = pf_stiff(E, A, Iz, L)
    return atba(r, k)


# Assemble structure stiffness matrix, contribution from one member. Related to structure
def pf_assemssm(imem, xy, conn, mprop, lm, ssm):
    """
    Superpose stiffness matrix of plane frame member on the structure
    stiffness matrix
    """
    K = pf_gstiff(imem, xy, conn, mprop)
    memdof = pf_get_dof(imem, conn, lm)
    for i in range(len(memdof)):
        if memdof[i]:
            for j in range(len(memdof)):
                if memdof[j]:
                    ii = memdof[i] - 1
                    jj = memdof[j] - 1
                    ssm[ii, jj] += K[i, j]
    return ssm


# Assemble structure stiffness matrix, contribution from all members. Related to structure
def pf_ssm(xy, conn, bc, mprop):
    """
    Assemble structure stiffness matrix by superposing stiffness matrices of
    individual members
    """
    n = xy.shape[0]
    lm, ndof = pf_calclm(n, bc)
    nmem = conn.shape[0]
    ssm = np.zeros((ndof, ndof), dtype=float)
    for imem in range(1, nmem + 1):
        ssm = pf_assemssm(imem, xy, conn, mprop, lm, ssm)
    return ssm, lm, ndof


# Assemble structure loadvector due to all joint loads. Related to structure
def pf_assem_loadvec_jl(lm, jtloads, P):
    """
    Superpose joint loads on structure load vector
    """
    nloads = jtloads.shape[0]
    for iload in range(nloads):
        jt = int(jtloads[iload, 0])
        jtdof = lm[jt - 1, :]
        for j in range(3):
            if jtdof[j]:
                i = jtdof[j] - 1
                P[i] += jtloads[iload, j + 1]
    return P


# Assemble structure loadvector due to one member load. Related to structure
def pf_assem_loadvec_ml(iload, xy, conn, lm, memloads, P):
    """
    Superpose equivalent joint loads due to loads applied on members onto
    structure load vector
    """
    imem = int(memloads[iload - 1, 0])
    xy1, xy2 = pf_get_endcoord(imem, xy, conn)
    L, dc = pf_calclen(xy1, xy2)
    r = pf_calcrot(dc)
    ml = memloads[iload - 1, 1:7]
    ml = ml.reshape(len(ml), 1)
    am = np.dot(-r.T, ml)
    memdof = pf_get_dof(imem, conn, lm)
    for i in range(6):
        if memdof[i]:
            ii = memdof[i] - 1
            P[ii] += am[i]
    return P


# Assemble structure loadvector due to all loads. Related to structure
def pf_loadvec(xy, conn, jtloads, memloads, ndof, lm):
    """
    Assemble structure load vector due to joint loads and loads applied
    directly on members
    """
    P = np.zeros((ndof, 1), dtype=float)
    P = pf_assem_loadvec_jl(lm, jtloads, P)
    nml = memloads.shape[0]
    for iload in range(1, nml + 1):
        P = pf_assem_loadvec_ml(iload, xy, conn, lm, memloads, P)
    return P


# Element end forces in one element. Related to element
def pf_mem_endforces(imem, xy, conn, mprop, memloads, lm, x):
    """
    Calculate member end forces from joint displacements, for one chosen member
    """
    xy1, xy2 = pf_get_endcoord(imem, xy, conn)
    L, dc = pf_calclen(xy1, xy2)
    E, A, Iz, L, dc = pf_get_memprop(imem, xy, conn, mprop)
    r = pf_calcrot(dc)
    u = np.zeros((6, 1), dtype=float)
    memdof = pf_get_dof(imem, conn, lm)
    for i in range(6):
        if memdof[i]:
            idof = memdof[i]
            u[i] = x[idof - 1]
    uu = np.dot(r, u)
    k = pf_stiff(E, A, Iz, L)
    f = np.zeros((6, 1), dtype=float)
    f = np.dot(k, uu)

    nml = memloads.shape[0]
    for i in range(nml):
        if memloads[i, 0] == imem:
            f += memloads[i, 1:].reshape(6, 1)
    return f


def pf_print_disp(lm, x):
    n = lm.shape[0]
    for i in range(n):
        print(f"{(i + 1):5d}", end="")
        for j in range(3):
            dof = int(lm[i, j])
            if dof:
                print(f"{x[dof - 1, 0]:12.6f}", end="")
            else:
                print(f"{0:12.6f}", end="")
        print()
    return 0


def input_data():
    # Weaver & Gere
    xy = np.array([[100.0, 75.0], [0.0, 75.0], [200.0, 0.0]])
    conn = np.array([[2, 1, 1], [1, 3, 1]])
    bc: Matrix = np.array([[2, 1, 1, 1], [3, 1, 1, 1]])
    mprop: Matrix = np.array([[1.0e4, 10, 1.0e3]], dtype=np.float_)
    jtloads: Matrix = np.array([[1, 0, -10.0, -1000.0]], dtype=np.float_)
    memloads = np.array([
        [1, 0.0, 12.0,200.0, 0.0, 12.0, -200.0],
        [2, -6.0, 8.0, 250.0, -6.0, 8.0, -250.0]], dtype=np.float_)

    # Hall & Kabaila
    # xy = np.array([[0.0, 8.0], [8.0, 8.0], [0.0, 4.0], [8.0, 4.0], [0.0, 0.0], [8.0, 0.0]], dtype=np.float_)
    # conn = np.array([[5, 3, 1], [3, 1, 1], [6, 4, 1], [4, 2, 1], [3, 4, 2], [1, 2, 2]], dtype=np.float_)
    # bc = np.array([[5, 1, 1, 1], [6, 1, 1, 1]], dtype=np.float_)
    # mprop = np.array([[20.0e6, 0.1, 0.0008], [20.0e6, 0.15, 0.0030]], dtype=np.float_)
    # jtloads = np.array([], dtype=np.float_)
    # memloads:Matrix = np.array([
    #     [6, 0.0, 10.0, 20.0, 0.0, 10.0, -20.0],
    #     [5, 0.0, 40.0, 53.3, 0.0, 40.0, -53.3]
    # ], dtype=np.float_)

    print_mat("\nJoint Coordinates\n", xy)
    print_mat("\nMember Connectivity\n", conn)
    print_mat("\nZero Boundary Conditions\n", bc)
    print_mat("\nMember Properties\n", mprop)
    return xy, conn, bc, mprop, jtloads, memloads


def main(xy, conn, bc, mprop, jtloads, memloads):
    ssm, lm, ndof = pf_ssm(xy, conn, bc, mprop)
    print_mat(f"\nNumber of degrees of freedom: {ndof:d}\n", lm)
    print_mat("\nStructure Stiffness Matrix\n", ssm)

    P = pf_loadvec(xy, conn, jtloads, memloads, ndof, lm)
    print_mat("\nLoad Vector\n", P)

    x = solve(ssm, P, assume_a="pos")
    print("\nDisplacements")
    pf_print_disp(lm, x * 1.0e3)

    print("\nMember End Forces")
    for imem in range(1, conn.shape[0] + 1):
        f = pf_mem_endforces(imem, xy, conn, mprop, memloads, lm, x)
        print_mat("%5d" % imem, f.T)
    return 0


# Main function
if __name__ == "__main__":
    xy, conn, bc, mprop, jtloads, memloads = input_data()
    main(xy, conn, bc, mprop, jtloads, memloads)
    # cProfile.run('main(xy, conn, bc, mprop, jtloads, memloads)', 'pf_profile')
    # import pstats
    # p = pstats.Stats('pf_profile')
    # p.strip_dirs().sort_stats('calls')
    # p.print_stats()
