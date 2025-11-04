import meshio

temp = meshio.read("../NumericalCases/Projects/qizhongji/last/MQ1330_remesh_resort.vtu")
with open("../output/cells.txt", 'w') as f:
    for key, eles in temp.cells_dict.items():
        for ele in eles:
            line = ",".join([str(node) for node in ele])
            f.write(f"{line}\n")
