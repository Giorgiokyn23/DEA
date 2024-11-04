import numpy as np
import scipy.linalg
import copy
import tqdm
from DEA import Circuit
import DEA




class model_pqc(Circuit):  # Eredita dalla tua classe Circuit
    def __init__(self, n_qubits, init=None):
        super().__init__(n_qubits, init)  # Chiamata al costruttore della classe base

    def conv(self, qubit0, qubit1, params):
        """Implement the conv operation using RX, RZ, and CNOT gates."""
        self.apply(self.RZ(qubit1))  # RZ(params[0]) on qubit1
        self.apply(self.CNOT(qubit0, qubit1))  # CNOT between qubit0 and qubit1
        self.apply(self.RZ(qubit0))  # RZ(params[1]) on qubit0
        self.apply(self.RX(qubit0))  # RX(params[2]) on qubit0
        self.apply(self.RZ(qubit0))  # RZ(params[3]) on qubit0
        self.apply(self.RZ(qubit1))  # RZ(params[4]) on qubit1
        self.apply(self.RX(qubit1))  # RX(params[5]) on qubit1
        self.apply(self.RZ(qubit1))  # RZ(params[6]) on qubit1
        self.apply(self.CNOT(qubit0, qubit1))  # CNOT between qubit0 and qubit1
        self.apply(self.RX(qubit0))  # RX(params[7]) on qubit0
        self.apply(self.RZ(qubit1))  # RZ(params[8]) on qubit1
        self.apply(self.CNOT(qubit0, qubit1))  # CNOT between qubit0 and qubit1
        self.apply(self.RZ(qubit0))  # RZ(params[9]) on qubit0
        self.apply(self.RX(qubit0))  # RX(params[10]) on qubit0
        self.apply(self.RZ(qubit0))  # RZ(params[11]) on qubit0
        self.apply(self.RZ(qubit1))  # RZ(params[12]) on qubit1
        self.apply(self.RX(qubit1))  # RX(params[13]) on qubit1
        self.apply(self.RZ(qubit1))  # RZ(params[14]) on qubit1

    def pool(self, sinkbit, targetbit, params):
        """Implement the pooling operation with rotations and CNOT."""
        self.apply(self.RZ(sinkbit))  # RZ(params[0]) on sinkbit
        self.apply(self.RY(sinkbit))  # RY(params[1]) on sinkbit
        self.apply(self.RZ(sinkbit))  # RZ(params[2]) on sinkbit
        self.apply(self.RZ(targetbit))  # RZ(params[3]) on targetbit
        self.apply(self.RY(targetbit))  # RY(params[4]) on targetbit
        self.apply(self.RZ(targetbit))  # RZ(params[5]) on targetbit
        self.apply(self.CNOT(sinkbit, targetbit))  # CNOT between sinkbit and targetbit
        self.apply(self.RZ(targetbit))  # RZ(params[6]) on targetbit
        self.apply(self.RY(targetbit))  # RY(params[7]) on targetbit
        self.apply(self.RZ(targetbit))  # RZ(params[8]) on targetbit

    def layer(self, params, wires):
        """Implement a layer with the conv and pool operations."""
        for idx, wire in enumerate(wires):
            if idx % 2 == 0 and idx < len(wires) - 1:
                self.conv(wire, wire + 1, params[:15])
        for idx, wire in enumerate(wires):
            if idx % 2 == 1 and idx < len(wires) - 1:
                self.conv(wire, wire + 1, params[:15])
        for idx, wire in enumerate(wires[:len(wires) // 2]):
            self.pool(wire, wire + len(wires) // 2, params[15:24])

    def layer2(self, params, wires):
        """Implement the second layer with conv and pool operations."""
        for idx, wire in enumerate(wires):
            if idx % 2 == 0 and idx < len(wires) - 1:
                self.conv(wire, wire + 1, params[:15])
        for idx, wire in enumerate(wires):
            if idx % 2 == 1 and idx < len(wires) - 1:
                self.conv(wire, wire + 1, params[:15])
        for idx, wire in enumerate(wires[:len(wires) // 2]):
            self.pool(wire, wire + len(wires) // 2, params[15:24])

    def build(self, params):
        """Build the entire PQC model."""
        self.layer(params[:24], range(self.n_qubits))  # Apply the first layer
        self.layer2(params[24:48], range(self.n_qubits // 2))  # Apply the second layer

# Esempio di utilizzo
n_qubits = 4
params = np.random.random(48)  # Inizializzazione dei parametri

C = model_pqc(n_qubits)
C.build(params)  # Costruisce il circuito e traccia i parametri


#Esegui IDEA sul circuito PQC
phi_dict = {
    0: [0, 15, 30],    # Parametri e gate associati
    1: [1, 16, 31],
    2: [2, 17, 32],
    3: [3, 18, 33],
    4: [4, 19, 34],
    5: [5, 20, 35],
    6: [6, 21, 36],
    7: [7, 22, 37],
    8: [8, 23, 38],
    9: [9, 24, 39],
    10: [10, 25, 40],
    11: [11, 26, 41],
    12: [12, 27, 42],
    13: [13, 28, 43],
    14: [14, 29, 44],
    15: [45,54],
    16: [46,55],
    17: [47,56],
    18: [48,57],
    19: [49,58],
    20: [50,59],
    21: [51,60],
    22: [52,61],
    23: [53,62],
    24: [63],
    25: [64],
    26: [65],
    27: [66],
    28: [67],
    29: [68],
    30: [69],
    31: [70],
    32: [71],
    33: [72],
    34: [73],
    35: [74],
    36: [75],
    37: [76],
    38: [77],
    39: [78],
    40: [79],
    41: [80],
    42: [81],
    43: [82],
    44: [83],
    45: [84],
    46: [85],
    47: [86],
}


A = DEA.IDEA(C,phi_dict, n_points=1)
print(A)
'''
if unique1:
    print("The list of independent parameters is:", indep_params1)
else:
    print("False:", indep_params1)
'''
