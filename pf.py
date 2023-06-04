import sqlite3
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
import tomllib
from scipy.linalg import solve

NDArrayInt = npt.NDArray[np.int_]
NDArrayFloat = npt.NDArray[np.float_]
Matrix = npt.NDArray[np.float_]


def atba(a, b):
    return np.dot(np.dot(a.T, b), a)


def array2df(x, columns=None, dtypes=None):
    df = pd.DataFrame(x, columns=columns)
    if dtypes:
        df = df.astype(dtypes)
    return df


# Location Matrix. Related to structure
def pf_calclm(n, df_bc: pd.DataFrame):
    """
    Calculate location matrix containing degree of freedom numbers of each node
    """
    lm = np.zeros((n, 3), dtype=np.int_)
    nd = 0
    ns = len(bc)
    for i in range(ns):
        node = int(df_bc.iloc[i, 0])  # type: ignore
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
def pf_stiff(E: float, A: float, Iz: float, L: float) -> NDArrayFloat:
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
    Calculate rotation matrix of a plane frame member given its direction cosines
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
def pf_get_endjts(imem, df_conn: pd.DataFrame):
    """Return numbers of first and second end of a plane frame member"""
    jt1 = df_conn.iloc[imem - 1, 0]
    jt2 = df_conn.iloc[imem - 1, 1]
    return jt1, jt2


# Utility function to get coordinates of end joints. Related to element
def pf_get_endcoord(imem, df_xy, df_conn):
    """Return the coordinates of the first and second end of a plane frame member"""
    jt1, jt2 = pf_get_endjts(imem, df_conn)
    xy1 = df_xy.iloc[int(jt1) - 1, :]  # type: ignore
    xy2 = df_xy.iloc[int(jt2) - 1, :]  # type: ignore
    return xy1, xy2


# Utility function to get member properties. Related to element
def pf_get_memprop(imem, df_xy, df_conn, df_mprop):
    """Return the properties of specified plane frame member"""
    m = df_conn.iloc[imem - 1, 2] - 1
    E = df_mprop.iloc[m, 0]
    A = df_mprop.iloc[m, 1]
    Iz = df_mprop.iloc[m, 2]
    xy1, xy2 = pf_get_endcoord(imem, df_xy, df_conn)
    L, dc = pf_calclen(xy1, xy2)
    return E, A, Iz, L, dc


# Utility function to get dof numbers at end joints. Related to element
def pf_get_dof(imem, df_conn, lm):
    """Return degree of freedom numbers of the ends of a plane frame member"""
    jt1, jt2 = pf_get_endjts(imem, df_conn)
    memdof = np.array([0, 0, 0, 0, 0, 0])
    memdof[0:3] = lm[jt1 - 1, :]  # type: ignore
    memdof[3:6] = lm[jt2 - 1, :]  # type: ignore
    return memdof


# Stiffness matrix in structure coordinate system. Related to element
def pf_gstiff(imem, df_xy, df_conn, df_mprop):
    """
    Calculate the stiffness matrix of a plane frame member in strucutre
    coordinate system
    """
    E, A, Iz, L, dc = pf_get_memprop(imem, df_xy, df_conn, df_mprop)
    r = pf_calcrot(dc)
    k = pf_stiff(E, A, Iz, L)
    return atba(r, k)


# Assemble structure stiffness matrix, contribution from one member. Related to structure
def pf_assemssm(imem, df_xy, df_conn, df_mprop, lm, ssm):
    """
    Superpose stiffness matrix of plane frame member on the structure
    stiffness matrix
    """
    K = pf_gstiff(imem, df_xy, df_conn, df_mprop)
    memdof = pf_get_dof(imem, df_conn, lm)
    for i in range(len(memdof)):
        if memdof[i]:
            for j in range(len(memdof)):
                if memdof[j]:
                    ii = memdof[i] - 1
                    jj = memdof[j] - 1
                    ssm[ii, jj] += K[i, j]
    return ssm


# Assemble structure stiffness matrix, contribution from all members. Related to structure
def pf_ssm(df_xy, df_conn, df_bc, df_mprop):
    """
    Assemble structure stiffness matrix by superposing stiffness matrices of
    individual members
    """
    n = len(df_xy)
    lm, ndof = pf_calclm(n, df_bc)
    nmem = len(conn)
    ssm = np.zeros((ndof, ndof), dtype=float)
    for imem in range(1, nmem + 1):
        ssm = pf_assemssm(imem, df_xy, df_conn, df_mprop, lm, ssm)
    return ssm, lm, ndof


# Assemble structure loadvector due to all joint loads. Related to structure
def pf_assem_loadvec_jl(lm, df_jtloads, P):
    """
    Superpose joint loads on structure load vector
    """
    nloads = len(jtloads)
    for iload in range(nloads):
        jt = int(df_jtloads.iloc[iload, 0])
        jtdof = lm[jt - 1, :]
        for j in range(3):
            if jtdof[j]:
                i = jtdof[j] - 1
                P[i] += df_jtloads.iloc[iload, j + 1]
    return P


# Assemble structure loadvector due to one member load. Related to structure
def pf_assem_loadvec_ml(iload, df_xy, df_conn, lm, df_memloads, P):
    """
    Superpose equivalent joint loads due to loads applied on members onto
    structure load vector
    """
    imem = int(df_memloads.iloc[iload - 1, 0])
    xy1, xy2 = pf_get_endcoord(imem, df_xy, df_conn)
    L, dc = pf_calclen(xy1, xy2)
    r = pf_calcrot(dc)
    ml = memloads[iload - 1, 1:7]
    ml = ml.reshape(len(ml), 1)
    am = np.dot(-r.T, ml)
    memdof = pf_get_dof(imem, df_conn, lm)
    for i in range(6):
        if memdof[i]:
            ii = memdof[i] - 1
            P[ii] += am[i]
    return P


# Assemble structure loadvector due to all loads. Related to structure
def pf_loadvec(df_xy, df_conn, df_jtloads, df_memloads, ndof, lm):
    """
    Assemble structure load vector due to joint loads and loads applied
    directly on members
    """
    P = np.zeros((ndof, 1), dtype=float)
    P = pf_assem_loadvec_jl(lm, df_jtloads, P)
    nml = len(df_memloads)
    for iload in range(1, nml + 1):
        P = pf_assem_loadvec_ml(iload, df_xy, df_conn, lm, df_memloads, P)
    return P


# Element end forces in one element. Related to element
def pf_mem_endforces(imem, df_xy, df_conn, df_mprop, df_memloads, lm, x):
    """
    Calculate member end forces from joint displacements, for one chosen member
    """
    xy1, xy2 = pf_get_endcoord(imem, df_xy, df_conn)
    L, dc = pf_calclen(xy1, xy2)
    E, A, Iz, L, dc = pf_get_memprop(imem, df_xy, df_conn, df_mprop)
    r = pf_calcrot(dc)
    u = np.zeros((6, 1), dtype=float)
    memdof = pf_get_dof(imem, df_conn, lm)
    for i in range(6):
        if memdof[i]:
            idof = memdof[i]
            u[i] = x[idof - 1]
    uu = np.dot(r, u)
    k = pf_stiff(E, A, Iz, L)
    f = np.zeros((6, 1), dtype=float)
    f = np.dot(k, uu)

    nml = len(df_memloads)
    for i in range(nml):
        if df_memloads.iloc[i, 0] == imem:
            print("***", df_memloads.iloc[i, 1:])
            f += df_memloads.iloc[i, 1:].values.reshape(6, 1)
    return f


def pf_print_disp(lm, x):
    n = len(lm)
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
    conn = np.array([[2, 1, 1], [1, 3, 1]], dtype=np.int_)
    bc: NDArrayFloat = np.array([[2, 1, 1, 1], [3, 1, 1, 1]])
    mprop: NDArrayFloat = np.array([[1.0e4, 10, 1.0e3]], dtype=np.float_)
    jtloads: Matrix = np.array([[1, 0, -10.0, -1000.0]], dtype=np.float_)
    memloads = np.array(
        [[1, 0.0, 12.0, 200.0, 0.0, 12.0, -200.0], [2, -6.0, 8.0, 250.0, -6.0, 8.0, -250.0]], dtype=np.float_
    )

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


def main(title, df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads):
    print(f"{title}\n{'='*len(title)}")
    ssm, lm, ndof = pf_ssm(df_xy, df_conn, df_bc, df_mprop)
    print_mat(f"\nNumber of degrees of freedom: {ndof:d}\n", lm)
    print_mat("\nStructure Stiffness Matrix\n", ssm)

    P = pf_loadvec(df_xy, df_conn, df_jtloads, df_memloads, ndof, lm)
    print_mat("\nLoad Vector\n", P)

    x = solve(ssm, P, assume_a="pos")
    print("\nDisplacements")
    pf_print_disp(lm, x * 1.0e3)

    print("\nMember End Forces")
    for imem in range(1, len(df_conn) + 1):
        f = pf_mem_endforces(imem, df_xy, df_conn, df_mprop, df_memloads, lm, x)
        print_mat("%5d" % imem, f.T)
        disp = array2df(f.T, ["Px1", "Py1", "Mz1", "Px2", "Py2", "Mz2"])
        print(disp)
    return 0


def read_toml(fname: str):
    with open(fname) as f:
        s = f.read()
    data = tomllib.loads(s)
    title = data["title"]
    xy = np.array(data["coordinates"]["xy"], dtype=np.float_)
    conn = np.array(data["connectivity"]["conn"], dtype=np.int_)
    bc = np.array(data["boundary"]["bc"], dtype=np.int_)
    mprop = np.array(data["materials"]["mprop"], dtype=np.float_)
    jtloads = np.array(data["jointloads"]["jtloads"], dtype=np.float_)
    memloads = np.array(data["memberloads"]["memloads"], dtype=np.float_)
    return title, xy, conn, bc, mprop, jtloads, memloads


def data2df(
    xy: NDArrayFloat,
    conn: NDArrayInt,
    bc: NDArrayInt,
    mprop: NDArrayFloat,
    jtloads: NDArrayFloat,
    memloads: NDArrayFloat,
):  # , conn, bc, mprop, jtloads, memloads
    df_xy = pd.DataFrame(xy, columns=["x", "y"])
    df_xy = df_xy.astype(np.float_)
    df_conn = pd.DataFrame(conn, columns=["node1", "node2", "mprop"])
    df_conn = df_conn.astype(int)
    df_bc = pd.DataFrame(bc, columns=["node", "ux", "uy", "rz"])
    df_bc = df_bc.astype(np.int_)
    df_mprop = pd.DataFrame(mprop, columns=["E", "A", "Iz"])
    df_mprop = df_mprop.astype(np.float_)
    df_jtloads = pd.DataFrame(jtloads, columns=["node", "Px", "Py", "Mz"])
    df_jtloads = df_jtloads.astype({"node": np.int_})
    df_memloads = pd.DataFrame(memloads, columns=["member", "Px1", "Py1", "Mz1", "Px2", "Py2", "Mz2"])
    df_memloads = df_memloads.astype({"member": np.int_})
    return df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads


def sqlite2df(dbfile: str):
    p = Path(dbfile)
    print(p, p.exists())
    if not p.exists():
        raise FileNotFoundError
    try:
        con = sqlite3.connect(dbfile)
        xy = pd.read_sql("SELECT * FROM xy", con)
        xy.drop("id", axis=1, inplace=True)

        conn = pd.read_sql("SELECT * FROM conn", con)
        conn.drop("id", axis=1, inplace=True)

        bc = pd.read_sql("SELECT * FROM bc", con)
        bc.drop("id", axis=1, inplace=True)

        mprop = pd.read_sql("SELECT * FROM mprop", con)
        mprop.drop("id", axis=1, inplace=True)

        jtloads = pd.read_sql("SELECT * FROM jtloads", con)
        jtloads.drop("id", axis=1, inplace=True)

        memloads = pd.read_sql("SELECT * FROM memloads", con)
        memloads.drop("id", axis=1, inplace=True)

        return xy, conn, bc, mprop, jtloads, memloads
    except Exception as e:
        raise e


# Utility functions. Not directly related to DSM
def print_mat(header: str, k: npt.NDArray[Any]) -> None:
    """
    Print the array k, preceded by a header string. Each element of k is
    printed using the format string fmt
    """
    m, n = k.shape
    print(f"{header} Size: {m} x {n}")
    for i in range(m):
        for j in range(n):
            print(f"{k[i, j]:12.4f}", end="")
        print()


def print_df(header: str, df: pd.DataFrame) -> None:
    print(f"{header}\n{'='*len(header)}")
    print(df.to_string(index=False))
    print()


# Main function
if __name__ == "__main__":
    # xy, conn, bc, mprop, jtloads, memloads = input_data()
    title, xy, conn, bc, mprop, jtloads, memloads = read_toml("weaver.toml")
    # df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = data2df(xy, conn, bc, mprop, jtloads, memloads)
    df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = sqlite2df("weaver.sqlite3")
    main(title, df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads)
