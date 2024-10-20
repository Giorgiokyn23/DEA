import numpy as np
from DEA import Circuit, IDEA  # Import the Circuit class from the DEA module

class CHE_Circuit(Circuit):
    def __init__(self, n_qubits, init=None):
        super().__init__(n_qubits, init)  # Call the constructor of the base class

    def Hadamard(self, k):
        """Defines the non-parametric Hadamard gate for qubit k."""
        H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)  # Hadamard matrix
        # Apply the non-parametric Hadamard gate to the k-th qubit
        return self._apply_single_qubit_gate(H, k)

    def _apply_single_qubit_gate(self, gate, k):
        """Applies a 2x2 gate matrix to qubit k in the circuit."""
        if k == 0:
            result = gate
        else:
            result = self._ID()  # Identity for the qubits not involved
        for j in range(1, self.n_qubits):
            if j == k:
                result = self._tensor(gate, result)
            else:
                result = self._tensor(self._ID(), result)
        return result

    def apply_Hadamard_RZ(self, data):
        """Applies Hadamard and RZ to each qubit."""
        for idx, x in enumerate(data):
            #self.apply(self.Hadamard(idx))  # Apply Hadamard to each qubit
            #self.apply(self.RX(idx))
            self.apply(self.RZ(idx))  # Apply RZ to each qubit with value from data

    def apply_entangling_layer(self, data):
        """Applies CNOT and RZ gates based on the data."""
        x0, x1, x2, x3 = data[0], data[1], data[2], data[3]

        # Apply CNOT and RZ operations based on the data
        for idx in range(self.n_qubits - 1):
            if idx == (self.n_qubits - 4):  # idx=0
                self.apply(self.CNOT(idx, idx + 1))
                self.apply(self.RZ(idx + 1))  # Parametric RZ gate
                self.apply(self.CNOT(idx, idx + 1))

            if idx == (self.n_qubits - 3):  # idx=1
                self.apply(self.CNOT(idx, idx + 1))
                self.apply(self.RZ(idx + 1))  # Parametric RZ gate
                self.apply(self.CNOT(idx, idx + 1))

            if idx == (self.n_qubits - 2):  # idx=2
                self.apply(self.CNOT(idx, idx + 1))
                self.apply(self.RZ(idx + 1))  # Parametric RZ gate
                self.apply(self.CNOT(idx, idx + 1))
                self.apply(self.CNOT(idx - 2, idx + 1))
                self.apply(self.RZ(idx + 1))  # Parametric RZ gate
                self.apply(self.CNOT(idx - 2, idx + 1))
                self.apply(self.CNOT(idx - 1, idx + 1))
                self.apply(self.RZ(idx + 1))  # Parametric RZ gate
                self.apply(self.CNOT(idx - 1, idx + 1))

    def build(self, params):
        """Builds the entire CHE circuit using random parameters."""
        self.apply_Hadamard_RZ(params[:self.n_qubits])  # Apply Hadamard and RZ with the first parameters
        self.apply_entangling_layer(params[:self.n_qubits])  # Apply CNOT and RZ gates with the remaining parameters

# Example of using the CHE circuit and applying IDEA
n_qubits = 4  # Number of qubits
params = np.random.random(4)  # Initialize random parameters (4 total parameters for the circuit)

C = CHE_Circuit(n_qubits)
C.build(params)  # Build the circuit

# Initialize the phi_dict for this circuit
phi_dict = {
        0: [0, 4, 7],
        1: [1, 4, 5, 7, 8],
        2: [2, 5, 6],
        3: [3, 6, 7, 8]
}

# Run IDEA on the circuit
A = IDEA(C, phi_dict, n_points=10)
print(A)



