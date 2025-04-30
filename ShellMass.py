import numpy as np

# ==============================
# 1. 输入参数定义（可修改区域）
# ==============================
nodes = np.array([  # 四节点坐标 (x, y, z) [单位：米]
    [-1.0, -1.0, 0.0],  # 节点1
    [1.0, -1.0, 0.0],  # 节点2
    [1.0, 1.0, 0.0],  # 节点3
    [-1.0, 1.0, 0.0]  # 节点4
])
nodes = np.array([  # 四节点坐标 (x, y, z) [单位：米]
    [0, 0, 0.0],  # 节点1
    [1.0, 0, 0.0],  # 节点2
    [1.0, 1.0, 0.0],  # 节点3
    [0, 1.0, 0.0]  # 节点4
])
thickness = 0.1  # 壳厚度 [米]
density = 7850  # 材料密度 [kg/m³]
# density = 1  # 材料密度 [kg/m³]


# ==============================
# 2. 核心计算函数
# ==============================
def compute_shell_mass_matrix(nodes, thickness, density):
    # ------------------------------
    # 初始化参数
    # ------------------------------
    n_nodes = 4  # 四节点壳单元
    dof_per_node = 3  # 平动自由度 (ux, uy, uz)
    total_dof = n_nodes * dof_per_node
    mass_matrix = np.zeros((total_dof, total_dof))

    # ------------------------------
    # 定义2x2高斯积分点 (ξ, η, 权重)
    # ------------------------------
    gauss_points = [
        (-1 / np.sqrt(3), -1 / np.sqrt(3), 1.0),
        (1 / np.sqrt(3), -1 / np.sqrt(3), 1.0),
        (1 / np.sqrt(3), 1 / np.sqrt(3), 1.0),
        (-1 / np.sqrt(3), 1 / np.sqrt(3), 1.0)
    ]

    # ------------------------------
    # 高斯积分循环
    # ------------------------------
    for xi, eta, w in gauss_points:
        # 1. 计算形状函数及其导数
        N = 0.25 * np.array([
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta)
        ])

        dN_dxi = 0.25 * np.array([
            [eta - 1, xi - 1],
            [1 - eta, -xi - 1],
            [eta + 1, xi + 1],
            [-eta - 1, 1 - xi]
        ])

        # 2. 计算雅可比矩阵和行列式
        J = np.dot(dN_dxi.T, nodes[:, :2])  # 仅需x,y坐标计算面内雅可比
        det_J = np.abs(np.linalg.det(J))

        # 3. 构造形函数矩阵H (3x12)
        H = np.zeros((3, total_dof))
        for i in range(n_nodes):
            H[0, i * dof_per_node] = N[i]  # ux
            H[1, i * dof_per_node + 1] = N[i]  # uy
            H[2, i * dof_per_node + 2] = N[i]  # uz

        # 4. 计算积分点贡献
        rho_t = density * thickness  # 面密度
        dV = det_J * thickness * w  # 体积微元
        mass_contribution = H.T @ H * rho_t * dV

        # # 3. 构造12x1形函数向量 N_all
        # N_all = np.zeros((total_dof, 1))
        # for i in range(n_nodes):
        #     for j in range(dof_per_node):
        #         N_all[i * dof_per_node + j, 0] = N[i]
        #
        # # 4. 计算积分点贡献
        # rho_t = density * thickness  # 面密度
        # dV = det_J * thickness * w  # 体积微元
        # mass_contribution = N_all @ N_all.T * rho_t * dV

        # 5. 累加到总质量矩阵
        mass_matrix += mass_contribution

    return mass_matrix


# ==============================
# 3. 执行计算并验证
# ==============================
if __name__ == "__main__":
    M = compute_shell_mass_matrix(nodes, thickness, density)

    # 输出总质量验证
    total_mass = density * thickness * 4.0  # 面积=2x2=4 m²
    calculated_mass = np.sum(np.diag(M))
    # calculated_mass = np.sum(M)
    error = np.abs(calculated_mass - total_mass)

    # 打印结果
    print("一致质量矩阵 (12x12):\n", M)
    print("\n验证结果:")
    print(f"理论总质量: {total_mass:.2f} kg")
    print(f"计算总质量: {calculated_mass:.2f} kg")
    print(f"绝对误差: {error:.6f} kg (应接近0)")