from msa import pf

title, xy, conn, bc, mprop, jtloads, memloads = pf.read_toml("weaver.toml")
df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads = pf.sqlite2df("weaver.sqlite3")
pf.main(title, df_xy, df_conn, df_bc, df_mprop, df_jtloads, df_memloads)
