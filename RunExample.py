from main import *

if __name__ ==  "__main__":
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\tempdirectory\DanWeiQuaShell\DanWeiQuaShell.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\singleQuaShell.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\catibeam.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\banshou.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\BiaoZhunSuanli.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\solid\bridge.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\TestSparse.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\Circle12w.cdb"
    # input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\TestSparse.cdb"
    input_file = r"D:\WorkSpace\FEM\testcases\ANSYS\shell\TestSparse11wnode.cdb"
    MyPyFEM(pathlib.Path(input_file), check_model=False)
