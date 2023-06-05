import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import tomllib
from numpy.typing import NDArray

NDArrayInt = NDArray[np.int_]
NDArrayFloat = NDArray[np.float_]


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
    df_conn = df_conn.astype(np.int_)
    df_bc = pd.DataFrame(bc, columns=["node", "ux", "uy", "rz"])
    df_bc = df_bc.astype(np.int_)
    df_mprop = pd.DataFrame(mprop, columns=["E", "A", "Iz"])
    df_mprop = df_mprop.astype(np.float_)
    df_jtloads = pd.DataFrame(jtloads, columns=["node", "Px", "Py", "Mz"])
    df_jtloads = df_jtloads.astype({"node": np.int_})
    df_memloads = pd.DataFrame(memloads, columns=["node", "Px1", "Py1", "Mz1", "Px2", "Py2", "Mz2"])
    df_memloads = df_memloads.astype({"node": np.int_})
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


def array2df(x, columns=None, dtypes=None):
    df = pd.DataFrame(x, columns=columns)
    if dtypes:
        df = df.astype(dtypes)
    return df


if __name__ == "__main__":
    title, xy, conn, bc, mprop, jtloads, memloads = read_toml("weaver.toml")
    print(title)
    print("Node Coordinates\n", xy)
    print("Member Connectivity\n", conn)
    print("Boundary Conditions\n", bc)
    print("Material Properties\n", mprop)
    print("Nodal Loads\n", jtloads)
    print(memloads)
    df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = data2df(xy, conn, bc, mprop, jtloads, memloads)
    print(df_xy, "\n", df_conn, "\n", df_bc, "\n", df_mprop, "\n", df_jtloads, "\n", df_memloads)
    xy, conn, bc, mprop, jtloads, memloads = sqlite2df("weaver.sqlite3")
    print(xy)
    print(conn)
    print(bc)
    print(mprop)
    print(jtloads)
    print(memloads)
    print(array2df(xy, ["x", "y"], {"x": float, "y": float}))
    print(array2df(jtloads, ["node", "Px", "Py", "Mz"], {"node": int, "Px": float, "Py": float, "Mz": float}))
