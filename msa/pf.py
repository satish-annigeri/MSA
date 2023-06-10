"""Python package for analysis of skeletal structures using the direct stiffness method. The code is based
on the theory and flowcharts in Weaver, W. and Gere, J.M., *Matrix Analysis of Framed Structures*,
CBS Publishers, New Delhi, 1986. The package `pf` implements the analysis of plane frames.
"""

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


def array2df(x, columns=None, dtypes=None):
    df = pd.DataFrame(x, columns=columns)
    if dtypes:
        df = df.astype(dtypes)
    return df


# Location Matrix. Related to structure
def pf_calclm(n, df_bc: pd.DataFrame):
    """Calculate location matrix containing degree of freedom numbers of each node

    Parameters
    ----------
    n : int
        Number of nodes in the structure
    df_bc : pd.DataFrame
        Boundary condition for each degree of freedom of each node. 1 is restrained and 0 is unrestrained.

    Returns
    -------
    lm : NDArrayInt
         Integer array with degree of freedom number of each nodal displacement
    nd : int
         Number of degrees of freedom of the structure
    """
    lm = np.zeros((n, 3), dtype=np.int_)
    nd = 0
    ns = len(df_bc)
    for i in range(ns):
        node = int(df_bc.iloc[i, 0])  # type: ignore
        lm[node - 1, 0:3] = df_bc.iloc[i, 1:4]

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
    """Calculate stiffness of a plane frame member in local coordinate system

    Parameters
    ----------
    E : float
        Modulus of elasticity of material of the member
    A : float
        Cross-section area of the member
    Iz : float
        Second moment of area of the member about the axis of bending
    L : float
        Length of the member

    Returns
    -------
    NDArrayFloat
        Stiffness matrix of member about local axes
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
    """Calculate rotation matrix of a plane frame member given its direction cosines

    Parameters
    ----------
    dc : NDArrayFloat shape (2,)
        Direction cosines of the member with reference to x-y axes

    Returns
    -------
    NDArrayFloat shape (6, 6)
        Rotation matrix to transform local stiffness matrix to global stiffness matrix
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
    """Calculate length and direction cosines of a plane frame member

    Parameters
    ----------
    xy1 : NDArrayFloat shape (2,)
        (x, y) coordinates of node 1
    xy2 : NDArrayFloat shape (2,)
        (x, y) coordinates of node 2

    Returns
    -------
    float
        Length of the member
    """
    delta = xy2 - xy1
    L = np.sqrt(np.sum(delta**2))
    dc = delta / L
    return L, dc


# Utility function to get number of end joints. Related to element
def pf_get_endjts(imem, df_conn: pd.DataFrame):
    """Return numbers of first and second end of a plane frame member

    Parameters
    ----------
    imem : int
        Number of member
    df_conn : pd.DataFrame
        Connectivity matrix

    Returns
    -------
    jt1 : int
        Node1
    jt2 : int
        Node 2
    """
    jt1, jt2 = df_conn.iloc[imem - 1, 0:2]
    return jt1, jt2


# Utility function to get coordinates of end joints. Related to element
def pf_get_endcoord(imem, df_xy, df_conn):
    """Return the coordinates of the first and second end of a plane frame member

    Parameters
    ----------
    imem : int
        Number of the member
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members

    Returns
    -------
    xy1 : Float
        Coordinate of node1
    xy2 : Float
        Coordinate of node2
    """
    jt1, jt2 = pf_get_endjts(imem, df_conn)
    xy1 = df_xy.iloc[int(jt1) - 1, :]  # type: ignore
    xy2 = df_xy.iloc[int(jt2) - 1, :]  # type: ignore
    return xy1, xy2


# Utility function to get member properties. Related to element
def pf_get_memprop(imem, df_xy, df_conn, df_mprop):
    """Return the material properties of specified plane frame member

    Parameters
    ----------
    imem : int
        Number of the member
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    df_mprop : pd.DataFrame
        Material properties table

    Returns
    -------
    E : float
        Modulus of elasticity of the member
    A : float
        Cross-section area of member
    Iz : float
        Second moment of area of cross-section of the member about the axis of bending
    L : float
        Length of the member
    dc : NDArrayFloat
        Direction cosines of the member
    """
    m = df_conn.iloc[imem - 1, 2] - 1
    E = df_mprop.iloc[m, 0]
    A = df_mprop.iloc[m, 1]
    Iz = df_mprop.iloc[m, 2]
    xy1, xy2 = pf_get_endcoord(imem, df_xy, df_conn)
    L, dc = pf_calclen(xy1, xy2)
    return E, A, Iz, L, dc


# Utility function to get dof numbers at end joints. Related to element
def pf_get_dof(imem, df_conn, lm):
    """Return degree of freedom numbers of the ends of a plane frame member

    Parameters
    ----------
    imem : int
        Number of the member
    df_conn : pd.DataFrame
        Connectivity of members
    lm : NDArrayFloat
        Location matrix, specifying degree of freedom number for each node

    Returns
    -------
    NDArrayInt
        Degree of freedom numbers for node1 and node2
    """
    jt1, jt2 = pf_get_endjts(imem, df_conn)
    memdof = np.array([0, 0, 0, 0, 0, 0])
    memdof[0:3] = lm[jt1 - 1, :]
    memdof[3:6] = lm[jt2 - 1, :]
    return memdof


# Stiffness matrix in structure coordinate system. Related to element
def pf_gstiff(imem, df_xy, df_conn, df_mprop):
    """Calculate the stiffness matrix of a plane frame member in strucutre
    coordinate system

    Parameters
    ----------
    imem : int
        Number of the member
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    df_mprop : pd.DataFrame
        Material properties table

    Returns
    -------
    NDArray
        Global stiffness matrix of the member
    """
    E, A, Iz, L, dc = pf_get_memprop(imem, df_xy, df_conn, df_mprop)
    r = pf_calcrot(dc)
    k = pf_stiff(E, A, Iz, L)
    return r.T @ k @ r


# Assemble structure stiffness matrix, contribution from one member. Related to structure
def pf_assemssm(imem, df_xy, df_conn, df_mprop, lm, ssm):
    """Superpose stiffness matrix of plane frame member on the structure
    stiffness matrix

    Parameters
    ----------
    imem : int
        Number of the member
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    df_mprop : pd.DataFrame
        Material properties table
    lm : NDArrayFloat
        Location matrix, specifying degree of freedom number for each node
    ssm : NDArrayFloat
        Structure stiffness matrix, initially a zero matrix

    Returns
    -------
    NDArrayFloat
        Structure stiffness matrix, after superposition of global stiffness of all members
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
    """Assemble structure stiffness matrix by superposing stiffness matrices of
    individual members

    Parameters
    ----------
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    df_bc : pd.DataFrame
        Zero boundary conditions table, for nodes with one or more constrained degrees of freedom
    df_mprop : pd.DataFrame
        Material properties table

    Returns
    -------
    ssm : NDArrayFloat
        Structure stiffness matrix of the structure
    lm : NDArrayInt
        Location matrix, specifying degree of freedom number for each node
    ndof : int
        Number of degrees of freedom of the structure
    """
    n = len(df_xy)
    lm, ndof = pf_calclm(n, df_bc)
    nmem = len(df_conn)
    ssm = np.zeros((ndof, ndof), dtype=float)
    for imem in range(1, nmem + 1):
        ssm = pf_assemssm(imem, df_xy, df_conn, df_mprop, lm, ssm)
    return ssm, lm, ndof


# Assemble structure loadvector due to all joint loads. Related to structure
def pf_assem_loadvec_jl(lm, df_jtloads, P):
    """Superpose joint loads on structure load vector

    Parameters
    ----------
    lm : NDArrayInt
        Location matrix, specifying degree of freedom number for each node
    df_jtloads : pd.DataFrame
        Nodal loads table
    P : NDArrayFloat
        Load vector

    Returns
    -------
    NDArrayFloat
        Modified load vector after superposing all nodal loads
    """
    nloads = len(df_jtloads)
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
    """Superpose equivalent joint loads due to loads applied on members onto
    structure load vector

    Parameters
    ----------
    iload : int
        Number of member load
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    lm : NDArrayInt
        Location matrix, specifying degree of freedom number for each node
    df_memloads : pd.DataFrame
        Member load table
    P : Load vector
        _description_

    Returns
    -------
    NDArrayFloat
        Modified load vector after superposing member load number `iload`
    """
    imem = int(df_memloads.iloc[iload - 1, 0])
    xy1, xy2 = pf_get_endcoord(imem, df_xy, df_conn)
    L, dc = pf_calclen(xy1, xy2)
    r = pf_calcrot(dc)
    ml = df_memloads.iloc[iload - 1, 1:7]
    # ml = ml.reshape(len(ml), 1)
    # ml = ml.values
    # am = np.dot(-r.T, ml)
    am = -r.T @ ml
    memdof = pf_get_dof(imem, df_conn, lm)
    for i in range(6):
        if memdof[i]:
            ii = memdof[i] - 1
            P[ii] += am[i]
    return P


# Assemble structure loadvector due to all loads. Related to structure
def pf_loadvec(df_xy, df_conn, df_jtloads, df_memloads, ndof, lm):
    """Assemble structure load vector due to joint loads and loads applied
    directly on members

    Parameters
    ----------
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    df_jtloads : pd.DataFrame
        Nodal loads table
    df_memloads : pd.DataFrame
        Member load table
    ndof : int
        Number of degrees of freedom of the structure
    lm : NDArrayInt
        Location matrix, specifying degree of freedom number for each node

    Returns
    -------
    NDArrayFloat
        Modified load vector after superposing all member loads
    """
    P = np.zeros((ndof, 1), dtype=float)
    P = pf_assem_loadvec_jl(lm, df_jtloads, P)
    nml = len(df_memloads)
    for iload in range(1, nml + 1):
        P = pf_assem_loadvec_ml(iload, df_xy, df_conn, lm, df_memloads, P)
    return P


# Element end forces in one element. Related to element
def pf_mem_endforces(imem, df_xy, df_conn, df_mprop, df_memloads, lm, x):
    """Calculate member end forces from joint displacements, for one chosen member

    Parameters
    ----------
    imem : int
        Number of the member
    df_xy : pd.DataFrame
        Coordinates of nodes
    df_conn : pd.DataFrame
        Connectivity of members
    df_mprop : pd.DataFrame
        Material properties table
    df_memloads : pd.DataFrame
        Member load table
    lm : NDArrayInt
        Location matrix, specifying degree of freedom number for each node
    x : NDArrayFloat
        Displacement vector

    Returns
    -------
    NDArrayFloat
        Vector of member end forces
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
    # uu = np.dot(r, u)
    uu = r @ u
    k = pf_stiff(E, A, Iz, L)
    f = np.zeros((6, 1), dtype=float)
    f = k @ uu
    # f = np.dot(k, uu)

    nml = len(df_memloads)
    for i in range(nml):
        if df_memloads.iloc[i, 0] == imem:
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


def pf_print_memendforces(imem, f):
    print(f"{imem:6d}", end=" ")
    for j in range(len(f)):
        print(f"{f.T[0, j]:12.4f}", end=" ")
    print()


def main(df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads):
    ssm, lm, ndof = pf_ssm(df_xy, df_conn, df_bc, df_mprop)
    print_mat(f"\nNumber of degrees of freedom: {ndof:d}\n", lm)
    print_mat("\nStructure Stiffness Matrix\n", ssm)

    P = pf_loadvec(df_xy, df_conn, df_jtloads, df_memloads, ndof, lm)
    print_mat("\nLoad Vector\n", P)

    x = solve(ssm, P, assume_a="pos")
    print("\nNode Displacements")
    pf_print_disp(lm, x * 1.0e3)

    print("\nMember End Forces")
    for imem in range(1, len(df_conn) + 1):
        f = pf_mem_endforces(imem, df_xy, df_conn, df_mprop, df_memloads, lm, x)
        pf_print_memendforces(imem, f)


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
    if len(jtloads) > 0:
        df_jtloads = pd.DataFrame(jtloads, columns=["node", "Px", "Py", "Mz"])
        df_jtloads = df_jtloads.astype({"node": np.int_})
    else:
        df_jtloads = pd.DataFrame(jtloads)
    if len(memloads) > 0:
        df_memloads = pd.DataFrame(memloads, columns=["member", "Px1", "Py1", "Mz1", "Px2", "Py2", "Mz2"])
        df_memloads = df_memloads.astype({"member": np.int_})
    else:
        df_memloads = pd.DataFrame(memloads)
    return df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads


def sqlite2df(dbfile: str):
    p = Path(dbfile)
    if not p.exists():
        raise FileNotFoundError
    try:
        con = sqlite3.connect(dbfile)

        problem = pd.read_sql("SELECT * FROM problem", con)
        title = str(problem.loc[0, "title"])

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

        return title, xy, conn, bc, mprop, jtloads, memloads
    except Exception as e:
        raise e


# Utility functions. Not directly related to DSM
def print_title(title: str):
    print(f"{title}\n{'-'*len(title)}")


def print_mat(header: str, k: npt.NDArray[Any]) -> None:
    """Print the array k, preceded by a header string. Each element of k is
    printed using the format string fmt

    Parameters
    ----------
    header : str
        Header to be printed before printing the matrix
    k : npt.NDArray[Any]
        Matrix to be printed
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


def pf_check(xy, conn):
    # Check connectivity
    last_node = len(xy)
    for imem in range(len(conn)):
        n1, n2 = pf_get_endjts(imem, conn)
        if (n1 > last_node) or (n2 > last_node):
            raise ValueError


# Main function
if __name__ == "__main__":
    # title, xy, conn, bc, mprop, jtloads, memloads = read_toml("weaver.toml")
    # df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = data2df(xy, conn, bc, mprop, jtloads, memloads)
    # print_title(title)
    # main(df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads)

    title, df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = sqlite2df("weaver.sqlite3")
    print_title(title)
    main(df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads)
