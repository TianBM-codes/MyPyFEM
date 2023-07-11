import numpy as np
from pathlib import Path


def mtx2mat(mtx_path: str, debug: bool = True):
    """
    Converts the Abaqus .MTX file to 2D numpy.nd arrays

    Parameters
    ----------
    mtx_path : str / pathlib.Path         of Abaqus .mtx file.
    debug : whether output

    Returns
    -------
    matrix : numpy.ndarray
        2D matrix of the .mtx file.
    """

    mtx_data = np.genfromtxt(str(mtx_path), delimiter=",")
    dof_per_node_count = int(max(np.max(mtx_data[:, 1]), np.max(mtx_data[:, 3])))
    node_count = int(max(np.max(mtx_data[:, 0]), np.max(mtx_data[:, 2])))
    total_dof_count = node_count * dof_per_node_count

    if debug:
        print(f"MTX File Path: {mtx_path}")
        print(f"Number of DOF per Node = {dof_per_node_count}")
        print(f"Total number of nodes = {node_count}")
        print(f"Total number of DOFs = {total_dof_count}")

    matrix = np.zeros([total_dof_count, total_dof_count])
    for mtx_row in mtx_data:
        mat_row = int(dof_per_node_count * (mtx_row[0] - 1) + (mtx_row[1] - 1))
        mat_col = int(dof_per_node_count * (mtx_row[2] - 1) + (mtx_row[3] - 1))
        matrix[mat_row, mat_col] = mtx_row[4]
        matrix[mat_col, mat_row] = mtx_row[4]

    return matrix


# data_dir = "./../testcases/static/linear/solid/hexaC3D8/stiff/"
# stiff_mat = mtx2mat(data_dir+"Job-1_STIF1.mtx")
stiff_mat = mtx2mat("../testcases/static/linear/solid/hexaC3D8/cube/Job-1_STIF1.mtx")