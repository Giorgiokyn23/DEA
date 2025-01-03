import numpy as np
import scipy.linalg
import copy
import tqdm

class Circuit:
    def __init__ (self,n_qubits,init=None):
        self.n_qubits = n_qubits
        self._dim = 2**self.n_qubits

        self._gates = []
        self._param_to_gate = []

        set_default_init = True
        if init != None:
            set_default_init = False
            # testing initialization
            try:
                self._wavefunction = np.asmatrix(init).astype(np.complex128).reshape((self._dim,1))
            except:
                print(f"invalid initialization - cannot cast to ({self.dim},1) np.complex128 matrix")
                set_default_init = True

            norm = 0.
            for j in range(self._dim):
                norm += abs(self._wavefunction[j,0])**2
            if norm == 0:
                print("invalid initialization - initialized with zero vector")
                set_default_init = True
            else:
                self._wavefunction = self._wavefunction/np.sqrt(norm)

        if set_default_init:
            self._wavefunction = np.asmatrix(np.zeros((self._dim,1),dtype=np.complex128))
            self._wavefunction[0,0] = 1.

        self.__id = np.identity(2,dtype=np.complex128)
        self.__sx = np.matrix([[0,1],[1,0]])
        self.__sy = np.matrix([[0,-1j],[1j,0]])
        self.__sz = np.matrix([[1,0],[0,-1]])
        self.__h = np.matrix([[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]])
        return

    def reset(self):
        self._wave_function = np.asmatrix(np.zeros((self._dim,1),dtype=np.complex128))
        self._wavefunction[0,0] = 1.
        self._gates = []
        self._param_to_gate = []
        return

    def reset_random(self):
        self._wave_function = np.asmatrix(np.zeros((self._dim,1),dtype=np.complex128))
        isset = False
        while not isset:
            for j in range(self._dim):
                self._wavefunction[j,0] = np.random.normal()
            norm = 0.
            for j in range(self._dim):
                norm += abs(self._wavefunction[j,0])**2
            if norm >0:
                isset = True
        self._wavefunction = self._wavefunction/np.sqrt(norm)
        self._gates = []
        self._param_to_gate = []
        return
    
    def _tensor(self,A,B):
        a1,a2 = np.shape(A)
        b1,b2 = np.shape(B)
        c1 = a1*b1
        c2 = a2*b2
        C = np.matrix(np.zeros((c1,c2),dtype=np.complex128))
        for j in range(c1):
            for k in range(c2):
                jb = j%b1
                kb = k%b2
                ja = int((j-jb)/b1)
                ka = int((k-kb)/b2)
                C[j,k] = A[ja,ka] * B[jb,kb]
        return C

    def _ID(self,k=0):
        return np.identity(self._dim,dtype=np.complex128)
    
    def _ZERO(self,k=0):
        return np.matrix(np.zeros((self._dim,self._dim),dtype=np.complex128))
    
    def _sX(self,k):
        if k == 0:
            sx = self.__sx
        else:
            sx = self.__id
        for j in range(1,self.n_qubits):
            if j == k:
                sx = self._tensor(self.__sx,sx)
            else:
                sx = self._tensor(self.__id,sx)
        return sx

    def _sY(self,k):
        if k == 0:
            sy = self.__sy
        else:
            sy = self.__id
        for j in range(1,self.n_qubits):
            if j == k:
                sy = self._tensor(self.__sy,sy)
            else:
                sy = self._tensor(self.__id,sy)
        return sy

    def _sZ(self,k):
        if k == 0:
            sz = self.__sz
        else:
            sz = self.__id
        for j in range(1,self.n_qubits):
            if j == k:
                sz = self._tensor(self.__sz,sz)
            else:
                sz = self._tensor(self.__id,sz)
        return sz

    def _h(self,k):
        if k == 0:
            h = self.__h
        else:
            h = self.__id
        for j in range(1,self.n_qubits):
            if j == k:
                h = self._tensor(self.__h,h)
            else:
                h = self._tensor(self.__id,h)
        return h

    def RX(self,k):
        r = lambda angle: np.cos(angle/2)*self._ID()-1j*np.sin(angle/2)*self._sX(k)
        dr = lambda angle: np.cos(angle/2)*self._sX(k)-1j*np.sin(angle/2)*self._ID()
        return (r,dr)

    def RY(self,k):
        r = lambda angle: np.cos(angle/2)*self._ID()-1j*np.sin(angle/2)*self._sY(k)
        dr = lambda angle: np.cos(angle/2)*self._sY(k)-1j*np.sin(angle/2)*self._ID()
        return (r,dr)

    def RZ(self,k):
        r = lambda angle: np.cos(angle/2)*self._ID()-1j*np.sin(angle/2)*self._sZ(k)
        dr = lambda angle: np.cos(angle/2)*self._sZ(k)-1j*np.sin(angle/2)*self._ID()
        return (r,dr)

    def RG(self,gate):
        r = lambda angle: scipy.linalg.expm(-1j*angle/2 * gate)
        dr = lambda angle: gate * scipy.linalg.expm(-1j*angle/2 * gate)
        return (r,dr)

    def H(self,k):
        return self._h(k)

    def CNOT(self,c,t):
        return (self._ID()+self._sZ(c))/2 + (self._ID()-self._sZ(c))/2 * self._sX(t)

    def CG(self,gate,c,*t):
        return (self._ID()+self._sZ(c))/2 + (self._ID()-self._sZ(c))/2 * gate(t)

    def apply(self,*gate):
        self._gates.append(gate[0])
        if isinstance(gate[0],tuple):
            self._param_to_gate.append(len(self._gates)-1)
        return

    def CircuitFn(self,args):
        s = self._wavefunction
        for idx in range(len(self._gates)):
            g = self._gates[idx]
            if isinstance(g,tuple):
                s = g[0](args[self._param_to_gate.index(idx)]) * s
            else:
                s = g*s
        return s

    '''
    def DiffCircuitFn(self, j, args):
        # Inizializza la wavefunction
        s = self._wavefunction
        print(f"Initial wavefunction: {s}")
        print(f"Args: {args}")
        print(f"Param to gate map: {self._param_to_gate}")

        for idx in range(len(self._gates)):
            g = self._gates[idx]
            print(f"\nProcessing gate at index {idx}: {g}")

            if isinstance(g, tuple):
                try:
                    arg_idx = self._param_to_gate.index(idx)
                    print(f"arg_idx for gate {idx}: {arg_idx}")

                    if arg_idx == j:
                        # Se arg_idx è uguale a j, applica il gate derivato
                        print(f"Applying derivative gate for param {arg_idx} on gate {idx}")
                        s = g[1](args[arg_idx]) * s
                    else:
                        # Altrimenti applica il gate normale
                        print(f"Applying gate for param {arg_idx} on gate {idx}")
                        s = g[0](args[arg_idx]) * s
                except IndexError as e:
                    print(f"Warning: param index {arg_idx} out of bounds. Skipping. Error: {e}")
                    continue
            else:
                print(f"Applying non-parametric gate at index {idx}")
                s = g * s

        print(f"Final wavefunction: {s}")
        return s

    '''
    def DiffCircuitFn(self, j ,args):
        s = self._wavefunction
        for idx in range(len(self._gates)):
            g = self._gates[idx]
            if isinstance(g,tuple):
                arg_idx = self._param_to_gate.index(idx)
                if arg_idx == j:
                    s = g[1](args[arg_idx]) * s
                else:
                    s = g[0](args[arg_idx]) * s
            else:
                s = g*s
        return s




def DEA(circuit, tol = 10**(-10), n_points = 1000,verbose=True):
    n_params = len(circuit._param_to_gate)
    np.random.seed(42)
    independent_at_point = []

    if verbose:
        print("Running DEA")
    for idx in tqdm.tqdm(range(n_points)):
        # choose random point
        theta = 2*np.pi*np.random.random(n_params)
        # compute S_N matrix
        S = np.asmatrix(np.zeros((n_params,n_params)))
        for j in range(n_params):
            for k in range(n_params):
                S[j,k] = (circuit.DiffCircuitFn(j,theta).H*circuit.DiffCircuitFn(k,theta)).real

        # find independent parameters at point
        indep_params = []
        for p in range(n_params):
            current_S = np.asmatrix(np.zeros((len(indep_params)+1,len(indep_params)+1)))
            current_params = copy.copy(indep_params)
            current_params.append(p)
            for j in range(len(current_params)):
                for k in range(len(current_params)):
                    current_S[j,k] = S[current_params[j],current_params[k]]

            if min(np.linalg.eigvals(current_S)) > tol:
                indep_params.append(p)
        # check whether this is a known set of independent parameters
        known = False
        for known_set_idx in range(len(independent_at_point)):
            known_set = independent_at_point[known_set_idx][0]
            same = (len(known_set) == len(indep_params))
            if same:
                for j in range(len(known_set)):
                    same = same and ( known_set[j] == indep_params[j] )
            known = known or same
            if same:
                known_idx = known_set_idx
        if known:
            independent_at_point[known_idx].append(theta)
        else:
            independent_at_point.append([indep_params,theta])

    if len(independent_at_point) == 1:
        return True, independent_at_point[0][0]
    return False, independent_at_point

def compute_S_matrix(circuit, n_params, theta):
    """Compute the S matrix based on the circuit's parameters."""
    S= np.asmatrix(np.zeros((n_params,n_params)))
    for j in range(n_params):
        for k in range(n_params):
            S[j,k]=(circuit.DiffCircuitFn(j,theta).H * circuit.DiffCircuitFn(k,theta)).real
        return S


def compute_phi_dict(circuit):
    """
    Computes the phi_dict and phi_values from the given circuit, correctly handling parameters
    that are reused across multiple gates.

    Parameters:
    - circuit: The quantum circuit object.

    Returns:
    - phi_dict: A dictionary mapping each parameter to the gates it affects (can map to multiple gates).
    - phi_values: A list of parameter values (initialized to None).
    """
    phi_dict = {}  # Dictionary for mapping parameters to gates
    phi_values = []  # List for storing parameter values

    param_to_gate_map = {}  # Temporary map from parameters to gate indices
    gate_counter = 0  # Keeps track of the gate index (ignores non-parametric gates)

    for gate_idx, gate in enumerate(circuit._gates):
        if isinstance(gate, tuple):  # Only consider parametric gates
            param_idx = circuit._param_to_gate.index(gate_idx)  # Get the parameter index
            print(circuit._gates)
            if param_idx not in param_to_gate_map:
                param_to_gate_map[param_idx] = []  # Initialize a list if it's not present

            param_to_gate_map[param_idx].append(gate_counter)  # Add this gate index to the parameter
            phi_values.append(None)  # Placeholder for the parameter value (could store actual values later)

            gate_counter += 1  # Only increment for parametric gates
        else:
            print(f"Skipping non-parametric gate at index {gate_idx}: {gate}")

    # Convert param_to_gate_map into phi_dict
    phi_dict = {k: v for k, v in param_to_gate_map.items()}

    return phi_dict, phi_values



# Attempt 26/09/2024

def compute_T_matrix(S, phi_dict, phi_values):
    """
    Compute the T matrix by following the generalized procedure.

    Parameters:
    - S: The base matrix representing gate interactions (numeric matrix).
    - phi_dict: A dictionary mapping each parameter phi to the list of gates they affect.
    - phi_values: A list of the phi parameter values (current values of the trainable parameters).

    Returns:
    - T: The resulting matrix T computed from S and the phi_dict mapping.
    """

    n_params = len(phi_dict)  # Number of parameters in phi_dict
    T = np.asmatrix(np.zeros((n_params, n_params)))  # Initialize T matrix with numerical entries

    #print(f"phi_dict: {phi_dict}")  # Debug: Print the phi_dict
    #print(f"phi_values:{phi_values}")
    # Loop over all parameter pairs (i, j) to compute T[i, j]
    for i in range(n_params):
        gates_i = phi_dict[i]  # Gates affected by parameter phi_i
        #print(f"gates_i for parameter {i}: {gates_i}")  # Debug: Print gates affected by parameter i

        for j in range(n_params):
            gates_j = phi_dict[j]  # Gates affected by parameter phi_j
            #print(f"gates_j for parameter {j}: {gates_j}")  # Debug: Print gates affected by parameter j

            # Initialize first and second terms for the scalar product
            first_term = []
            second_term = []

            # First term: Iterate over gates in gates_i
            for n in gates_i:
                if n >= S.shape[0]:
                    raise IndexError(f"Index {n} is out of bounds for S with size {S.shape[0]}")
                #print(f"n in gates_i: {n}")  # Debug: Print n

                # If gate n has more parameters, construct the product term
                factors_i = [phi_values[p] for p in range(n_params) if n in phi_dict[p] and p != i]
                if factors_i:
                    product_term_i = np.prod(factors_i)
                    first_term.append((product_term_i, n))  # Append the product and the gate index
                else:
                    first_term.append((1, n))  # Append 1 if no other parameters are present

            # Second term: Iterate over gates in gates_j
            for k in gates_j:
                if k >= S.shape[1]:
                    raise IndexError(f"Index {k} is out of bounds for S with size {S.shape[1]}")
                #print(f"k in gates_j: {k}")  # Debug: Print k

                factors_j = [phi_values[p] for p in range(n_params) if k in phi_dict[p] and p != j]
                if factors_j:
                    product_term_j = np.prod(factors_j)
                    second_term.append((product_term_j, k))  # Append the product and the gate index
                else:
                    second_term.append((1, k))  # Append 1 if no other parameters are present

            # Now, compute the scalar product between the two terms
            T_ij = 0
            for term_i, n in first_term:
                for term_j, k in second_term:
                    #print(f"Calculating T_ij for gates {n}, {k}")  # Debug: Print n, k before calculation
                    T_ij += term_i * term_j * S[n, k]

            # Store the result in T[i, j]
            T[i, j] = T_ij

    return T


def find_independent_parameters(T, n_params, tol, independent_at_point, theta):
    """
    Find the independent parameters at a given point using the T matrix.

    Parameters:
    - T: The matrix representing parameter interactions.
    - n_params: Number of parameters in the circuit.
    - tol: Tolerance for considering eigenvalues as zero.
    - independent_at_point: List of independent parameters found so far.
    - theta: The current parameter values.

    Returns:
    - independent_at_point: Updated list of independent parameters.
    """
    indep_params = []

    # Find independent parameters
    for p in range(n_params):
        current_T = np.asmatrix(np.zeros((len(indep_params) + 1, len(indep_params) + 1)))
        current_params = copy.copy(indep_params)
        current_params.append(p)

        for j in range(len(current_params)):
            for k in range(len(current_params)):
                current_T[j, k] = T[current_params[j], current_params[k]]

        #print("smallest EV:",min(np.linalg.eigvals(current_T)))
        # Check if the minimum eigenvalue is greater than the tolerance (indicating independence)
        if min(np.linalg.eigvals(current_T)) > tol:
            indep_params.append(p)

    # Check if this is a known set of independent parameters
    known = False
    for known_set_idx in range(len(independent_at_point)):
        known_set = independent_at_point[known_set_idx][0]
        same = (len(known_set) == len(indep_params))
        if same:
            for j in range(len(known_set)):
                same = same and (known_set[j] == indep_params[j])
        known = known or same
        if same:
            known_idx = known_set_idx

    # Update the list of independent parameters
    if known:
        independent_at_point[known_idx].append(theta)
    else:
        independent_at_point.append([indep_params, theta])

    return independent_at_point


#Fixed the dimension bug
def IDEA(circuit, phi_dict, tol=10 ** (-10), n_points=1000, verbose=True):
    """
    Perform Integrated Dimensional Expressivity Analysis (DEA) using the computed T matrix.

    Parameters:
    - circuit: The quantum circuit being analyzed.
    - phi_dict: The dictionary mapping each parameter to the gates it affects.
    - n_params: Number of parameters in the circuit.
    - tol: Tolerance for considering eigenvalues.
    - n_points: Number of random points to sample.
    - verbose: If True, print progress information.

    Returns:
    - A tuple with:
        - A boolean indicating if the parameters are independent.
        - The set of independent parameters if they exist.
    """
    # Set a fixed seed for reproducibility
    np.random.seed(42)
    n_params = len(circuit._param_to_gate)
    independent_at_point = []

    if verbose:
        print("Running IDEA with T matrix")

    for idx in tqdm.tqdm(range(n_points)):
        # Choose random parameter values for phi_values
        #print(n_params)
        theta = 2 * np.pi * np.random.random(n_params)
        #print(f"Iteration {idx}, Random theta: {theta}")

        # Compute the S matrix for the given circuit and phi_values
        S_matrix = np.asmatrix(np.zeros((n_params, n_params)))
        for j in range(n_params):
            for k in range(n_params):
                S_matrix[j, k] = (circuit.DiffCircuitFn(j, theta).H * circuit.DiffCircuitFn(k, theta)).real

        # Print the computed S_matrix for debugging purposes
        #print(f"S_matrix at iteration {idx}: \n{S_matrix}")

        # Compute the T matrix using the computed S matrix, phi_dict, and phi_values

        T = compute_T_matrix(S_matrix, phi_dict, theta)
        dimT= len(phi_dict)
        # Call the modularized function to find independent parameters
        independent_at_point = find_independent_parameters(T, dimT, tol, independent_at_point, theta)

    # Return the result: if only one set of independent parameters was found, return it
    if len(independent_at_point) == 1:
        return True, independent_at_point[0][0]
    return False, independent_at_point




