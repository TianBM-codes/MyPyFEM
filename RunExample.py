from main import *

if __name__ ==  "__main__":
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\tempdirectory\DanWeiQuaShell\DanWeiQuaShell.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\singleQuaShell.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\catibeam.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\banshou.cdb"
    input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\BiaoZhunSuanli.cdb"

    #  Unable to allocate 3.67 TiB for an array with shape (710604, 710604) and data type float64
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\Circle12w.cdb"
    MyPyFEM(pathlib.Path(input_file))
