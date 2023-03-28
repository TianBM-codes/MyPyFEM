import numpy as np
from element.ElementBase import *


class TetraElement(ElementBaseClass):
    """ Tetra Element class """

    def __init__(self, eid=None):
        super().__init__(eid)
        self._nodes_count = 2  # Each element has 2 nodes
        self._nodes = [None for _ in range(self._nodes_count)]
        self._vtp_type = "tetra"
        # self._ND = 6
        # self._LocationMatrix = np.zeros(self._ND, dtype=np.int)

    def GenerateLocationMatrix(self):
        """
        Generate location matrix: the global equation number that
        corresponding to each DOF of the element
        """
        i = 0
        for N in range(self._nodes_count):
            for D in range(3):
                self._location_matrix[i] = self._nodes[N].bcode[D]
                i += 1

    def SizeOfStiffnessMatrix(self):
        """
        Return the size of the element stiffness matrix
        (stored as an array column by column)
        For 2 node bar element, element stiffness is a 6x6 matrix,
        whose upper triangular part has 21 elements
        """
        return 21

    def ElementStiffness(self, stiffness):
        """
        Calculate element stiffness matrix
        Upper triangular matrix, stored as an array column by colum
        starting from the diagonal element
        """
        for i in range(self.SizeOfStiffnessMatrix()):
            stiffness[i] = 0.0

        # Calculate bar length
        # dx = x2-x1, dy = y2-y1, dz = z2-z1
        DX = np.zeros(3)
        for i in range(3):
            DX[i] = self._nodes[1].XYZ[i] - self._nodes[0].XYZ[i]

        # Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
        DX2 = np.zeros(6)
        DX2[0] = DX[0] * DX[0]
        DX2[1] = DX[1] * DX[1]
        DX2[2] = DX[2] * DX[2]
        DX2[3] = DX[0] * DX[1]
        DX2[4] = DX[1] * DX[2]
        DX2[5] = DX[0] * DX[2]

        L2 = DX2[0] + DX2[1] + DX2[2]
        L = np.sqrt(L2)

        # Calculate element stiffness matrix
        material = self._element_material

        k = material.E * material.Area / L / L2

        stiffness[0] = k * DX2[0]
        stiffness[1] = k * DX2[1]
        stiffness[2] = k * DX2[3]
        stiffness[3] = k * DX2[2]
        stiffness[4] = k * DX2[4]
        stiffness[5] = k * DX2[5]
        stiffness[6] = k * DX2[0]
        stiffness[7] = -k * DX2[5]
        stiffness[8] = -k * DX2[3]
        stiffness[9] = -k * DX2[0]
        stiffness[10] = k * DX2[1]
        stiffness[11] = k * DX2[3]
        stiffness[12] = -k * DX2[4]
        stiffness[13] = -k * DX2[1]
        stiffness[14] = -k * DX2[3]
        stiffness[15] = k * DX2[2]
        stiffness[16] = k * DX2[4]
        stiffness[17] = k * DX2[5]
        stiffness[18] = -k * DX2[2]
        stiffness[19] = -k * DX2[4]
        stiffness[20] = -k * DX2[5]

    def ElementStress(self, stress, displacement):
        """
        Calculate element stress
        """
        material = self._element_material

        DX = np.zeros(3)
        L2 = 0

        for i in range(3):
            DX[i] = self._nodes[1].XYZ[i] - self._nodes[0].XYZ[i]
            L2 += (DX[i] * DX[i])

        S = np.zeros(6)
        for i in range(3):
            S[i] = -DX[i] * material.E / L2
            S[i + 3] = -S[i]

        stress[0] = 0.0
        for i in range(6):
            if self._LocationMatrix[i]:
                stress[0] += (S[i] * displacement[self._LocationMatrix[i] - 1])
