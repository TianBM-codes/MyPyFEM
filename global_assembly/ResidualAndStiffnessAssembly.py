import numpy as np
from femdb.NLFEMDataBase import NLFEMDataBase
from femdb.Material import *
from femdb.Plasticity import PlasticDeformationState
from femdb.ElementGroup import ElementGroup
from utils.GlobalEnum import *
from scipy.sparse import coo_matrix, csr_matrix


def ResidualAndStiffnessAssembly(grp: ElementGroup):
    """
    Computes and assemble residual force vector and global tangent stiffness
    matrix except surface (line) element pressure contributions.
    """
    """
    Initialism variants
    """
    fem_db = NLFEMDataBase()
    global_k = fem_db.global_k
    right_hand_item = fem_db.right_hand_item
    MESH = fem_db.Mesh

    grp_ele_info = grp.element_info
    n_dofs_elem = grp_ele_info.n_dofs_elem
    n_nodes_elem = grp_ele_info.n_nodes_elem
    ngauss = grp_ele_info.ngauss
    ndim = GetDomainDimension()
    T_int = fem_db.right_hand_item.T_int

    right_hand_item.external_load = fem_db.SolveControl.xlmax * right_hand_item.nominal_external_load

    """
    Pre-allocation memory to indexi, indexj and data for sparse assembly of
    the stiffness matrix.
    """
    n_components_mean_dilatation = fem_db.Material.n_nearly_incompressible * np.square(n_dofs_elem)
    n_components_displacement_based = (MESH.nelem * np.square(n_dofs_elem) * ngauss
                                       + MESH.nelem * np.square(n_nodes_elem) * ndim * ngauss)
    n_components = n_components_mean_dilatation + n_components_displacement_based

    """
    Initialise counter for storing sparse information into 
    global tangent stiffness matrix.
    """
    global_k.counter = 0
    global_k.indexi = np.zeros((n_components, ), dtype=np.uint32)
    global_k.indexj = np.zeros((n_components, ), dtype=np.uint32)
    global_k.stiffness = np.zeros((n_components, ), dtype=float)

    """
    Main element loop
    """
    for _, grp in fem_db.ElementGroupHash.items():
        for ele_idx in range(len(grp.eles)):
            ele = grp.eles[ele_idx]
            node_ids = ele.GetNodeSearchIndex()
            mat_id = ele.mat_id
            ele_id = ele.search_ele_idx
            mat = fem_db.Material[mat_id]
            if isinstance(mat, HyperElasticPlasticInPrincipal):
                plasticity_element = PlasticDeformationState()
                plasticity_element.epbar = grp.global_plasticity.epbar[:, ele_id]
                plasticity_element.invCp = grp.global_plasticity.invCp[:, :, :, ele_id]
            else:
                raise NoImplSuchMaterial(mat.GetName())

            xlocal = fem_db.Geom.x[:, node_ids]
            x0local = fem_db.Geom.x0[:, node_ids]
            Ve = fem_db.Geom.V_ele[ele_idx]
            from element_calculation.ElementForceAndStiffness import ElementForceAndStiffness
            T_internal, PLAST_element = ElementForceAndStiffness(xlocal, x0local, mat_id, Ve, ele, grp, ele_idx)

            """
            Assemble element contribution into global internal force vector. 
            Utility function for the assembly of element force vectors into the 
            global force vector. 
            """
            global_dofs = MESH.dof_nodes[:, ele.search_node_ids]
            T_int[global_dofs.flatten()] = T_int[global_dofs.flatten()] + T_internal

            """
            Storage of updated value of the internal variables.
            """
            grp.global_plasticity.invCp[:, :, :, ele_idx] = PLAST_element.invCp
            grp.global_plasticity.epbar[:, ele_idx] = PLAST_element.epbar

    """
    Global tangent stiffness matrix sparse assembly except pressure contributions. 
    """
    S = coo_matrix((global_k.stiffness, (global_k.indexi, global_k.indexj)), dtype=float)
    global_k.stiffness = csr_matrix(S)

    """
    Compute global residual force vector except pressure contributions.
    """
    right_hand_item.residual = right_hand_item.T_int - right_hand_item.external_load
