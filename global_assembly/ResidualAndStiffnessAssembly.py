import numpy as np
from femdb.NLFEMDataBase import NLFEMDataBase
from femdb.NLDomain import NLDomain
from utils.GlobalEnum import *
from scipy.sparse import coo_matrix, csr_matrix


def ResidualAndStiffnessAssembly():
    """
    Computes and assemble residual force vector and global tangent stiffness
    matrix except surface (line) element pressure contributions.
    """
    """
    Initialism variants
    """
    fem_db = NLFEMDataBase()
    nl_domain = NLDomain()
    global_k = nl_domain.global_k
    aux_variant = nl_domain.aux_variant
    right_hand_item = nl_domain.right_hand_item
    plastics = nl_domain.plastics

    n_dofs_elem = aux_variant.n_dofs_elem
    ngauss = aux_variant.ngauss
    ndim = GlobalInfor[GlobalVariant.Dimension]

    right_hand_item.external_load = fem_db.SolveControl.xlmax * right_hand_item.nominal_external_load

    """
    Pre-allocation memory to indexi, indexj and data for sparse assembly of
    the stiffness matrix.
    """
    n_components_mean_dilatation = fem_db.Material.n_nearly_incompressible * np.square(n_dofs_elem, 2)
    n_components_displacement_based = (fem_db.Mesh.nelem * np.square(n_dofs_elem, 2)
                                       + fem_db.Mesh.nelem * np.square(n_dofs_elem, 2) * ndim * ngauss)
    n_components = n_components_mean_dilatation + n_components_displacement_based

    """
    Initialise counter for storing sparse information into 
    global tangent stiffness matrix.
    """
    global_k.counter = 1
    global_k.indexi = np.zeros((n_components, 1))
    global_k.indexj = np.zeros((n_components, 1))
    global_k.stiffness = np.zeros((n_components, 1))

    """
    Main element loop
    """
    for _, grp in fem_db.ElementGroupHash.items():
        for ele in grp.eles:
            node_ids = ele.GetNodes()
            mat_id = ele.mat_id
            mat = fem_db.Material.Mat[mat_id]
            ele_id = ele.id_key
            epbar = plastics.plastic_deformation_state.epbar[:, ele_id]
            invCp = plastics.plastic_deformation_state.invCp[:, :, :, ele_id]
            # TODO 这里需要节点排好，对应关系另外存储，不要每次都查询dict, 只有在读取输入文件和写结果的时候初始化各种dict
            xlocal = fem_db.Geom.x[node_ids, :]
            x0local = fem_db.Geom.x0[node_ids, :]
            Ve = fem_db.Geom.V_ele[node_ids, :]
            from element_calculation.ElementForceAndStiffness import ElementForceAndStiffness
            ElementForceAndStiffness(xlocal, x0local, mat_id, Ve, ele)

            """
            Storage of updated value of the internal variables. 
            """
            # TODO: 首先将节点排好序后，再储存 plasticity

    """
    Global tangent stiffness matrix sparse assembly except pressure contributions. 
    """
    S = coo_matrix(global_k.indexi, global_k.indexj, global_k.stiffness)
    global_k.stiffness = csr_matrix(S)
    """
    Compute global residual force vector except pressure contributions.
    """
    right_hand_item.residual = right_hand_item.T_int - right_hand_item.external_load
