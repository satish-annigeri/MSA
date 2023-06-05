import numpy as np

from msa import pf


def input_data(prob: str):
    match prob.lower():
        case "weaver":
            # Weaver & Gere
            title = "Plane Frame (Weaver and Gere)"
            xy = np.array([[100.0, 75.0], [0.0, 75.0], [200.0, 0.0]])
            conn = np.array([[2, 1, 1], [1, 3, 1]], dtype=np.int_)
            bc = np.array([[2, 1, 1, 1], [3, 1, 1, 1]], dtype=np.int_)
            mprop = np.array([[1.0e4, 10, 1.0e3]], dtype=np.float_)
            jtloads = np.array([[1, 0, -10.0, -1000.0]], dtype=np.float_)
            memloads = np.array(
                [[1, 0.0, 12.0, 200.0, 0.0, 12.0, -200.0], [2, -6.0, 8.0, 250.0, -6.0, 8.0, -250.0]], dtype=np.float_
            )
        case "hall":
            # Hall & Kabaila
            title = "Plane Frame (Hall and Kabaila)"
            xy = np.array([[0.0, 8.0], [8.0, 8.0], [0.0, 4.0], [8.0, 4.0], [0.0, 0.0], [8.0, 0.0]], dtype=np.float_)
            conn = np.array([[5, 3, 1], [3, 1, 1], [6, 4, 1], [4, 2, 1], [3, 4, 2], [1, 2, 2]], dtype=np.int_)
            bc = np.array([[5, 1, 1, 1], [6, 1, 1, 1]], dtype=np.int_)
            mprop = np.array([[20.0e6, 0.1, 0.0008], [20.0e6, 0.15, 0.0030]], dtype=np.float_)
            jtloads = np.array([], dtype=np.float_)
            memloads = np.array(
                [[6, 0.0, 10.0, 20.0, 0.0, 10.0, -20.0], [5, 0.0, 40.0, 53.3, 0.0, 40.0, -53.3]], dtype=np.float_
            )
        case _:
            # Weaver & Gere
            title = "Plane Frame (Weaver and Gere)"
            xy = np.array([[100.0, 75.0], [0.0, 75.0], [200.0, 0.0]])
            conn = np.array([[2, 1, 1], [1, 3, 1]], dtype=np.int_)
            bc = np.array([[2, 1, 1, 1], [3, 1, 1, 1]], dtype=np.int_)
            mprop = np.array([[1.0e4, 10, 1.0e3]], dtype=np.float_)
            jtloads = np.array([[1, 0, -10.0, -1000.0]], dtype=np.float_)
            memloads = np.array(
                [[1, 0.0, 12.0, 200.0, 0.0, 12.0, -200.0], [2, -6.0, 8.0, 250.0, -6.0, 8.0, -250.0]], dtype=np.float_
            )
    return title, xy, conn, bc, mprop, jtloads, memloads


def echo_input(title, xy, conn, bc, mprop):
    pf.print_title(title)
    pf.print_mat("\nJoint Coordinates\n", xy)
    pf.print_mat("\nMember Connectivity\n", conn)
    pf.print_mat("\nZero Boundary Conditions\n", bc)
    pf.print_mat("\nMember Properties\n", mprop)


if __name__ == "__main__":
    title, xy, conn, bc, mprop, jtloads, memloads = input_data("weaver")
    pf.print_title(title)
    df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = pf.data2df(xy, conn, bc, mprop, jtloads, memloads)
    pf.main(df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads)

    print(f"\n{'='*90}\n")

    title, xy, conn, bc, mprop, jtloads, memloads = input_data("hall")
    pf.print_title(title)
    df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = pf.data2df(xy, conn, bc, mprop, jtloads, memloads)
    pf.main(df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads)

    # title, xy, conn, bc, mprop, jtloads, memloads = pf.read_toml("weaver.toml")
    # df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = pf.sqlite2df("weaver.sqlite3")
