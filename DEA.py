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
                self._wavefunction = np.asmatrix(init).astype(np.complex128).reshape((self.dim,1))
            except:
                print(f"invalid initialization - cannot cast to ({self.dim},1) np.complex128 matrix")
                set_default_init = True

            norm = 0.
            for j in range(self.dim):
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
        return

    def reset(self):
        self._wave_function = np.asmatrix(np.zeros((self.dim,1),dtype=np.complex128))
        self._wavefunction[0,0] = 1.
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

    def DiffCircuitFn(self,j,args):
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
            current_S = np.asmatrix(np.zeros((p+1,p+1)))
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
    Compute the phi_dict and phi_values for a given circuit object.

    Parameters:
    - circuit: The Circuit object that contains gates and parameter information.

    Returns:
    - phi_dict: A dictionary mapping each parameter index to the list of gate indices it affects.
    - phi_values: A list of the current phi parameter values (trainable parameters).
    """
    phi_dict = {}
    phi_values = []

    param_idx = 0  # Index for phi parameters
    for idx, gate in enumerate(circuit._gates):
        # Check if the gate is parametric by checking if it's a tuple (parameterized gate)
        if isinstance(gate, tuple):
            if param_idx not in phi_dict:
                phi_dict[param_idx] = []
            # Append the current gate index to the list of gates associated with the parameter
            phi_dict[param_idx].append(idx)
            # Store the parameter value for this gate (assuming each gate has its angle as a parameter)
            phi_values.append(None)  # Placeholder for now, as no specific values provided
            param_idx += 1

    return phi_dict, phi_values


# Attempt 26/09/2024
def compute_T_matrix(S, phi_dict, phi_values):
    """
    Compute the T matrix by following the generalized procedure.

    Parameters:
    - S: The base matrix representing gate interactions
    - phi_dict: A dictionary mapping each parameter phi to the list of gates they affect.
    - phi_values: A list of the phi parameter values (current values of the trainable parameters).

    Returns:
    - T: The resulting matrix T computed from S and the phi_dict mapping
    """

    n_params = len(phi_dict)  # Number of parameters in phi_dict

    # Initialize T matrix as an empty matrix using np.asmatrix to represent the matrix structure
    T = np.asmatrix(np.zeros((n_params, n_params)))

    # Loop over all parameter pairs (i, j) to compute T[i, j]
    for i in range(n_params):
        gates_i = phi_dict[i]  # Gates affected by parameter phi_i

        for j in range(n_params):
            gates_j = phi_dict[j]  # Gates affected by parameter phi_j

            # Initialize first and second terms for the scalar product
            first_term = []  # List to store terms related to parameter phi_i
            second_term = []  # List to store terms related to parameter phi_j

            # First term: Iterate over gates in gates_i
            for n in gates_i:
                # Check if gate n is shared with other parameters, i.e., apply the product rule
                factors_i = [phi_values[p] for p in range(n_params) if n in phi_dict[p] and p != i]

                if factors_i:
                    # If there are multiple parameters in the gate, compute their product
                    product_term_i = np.prod(factors_i)
                    first_term.append((product_term_i, n))  # Append (product, gate index n)
                else:
                    # If no other parameters, append 1 (since there's no additional multiplication)
                    first_term.append((1, n))

            # Second term: Iterate over gates in gates_j
            for k in gates_j:
                # Check if gate k is shared with other parameters, apply the product rule if needed
                factors_j = [phi_values[p] for p in range(n_params) if k in phi_dict[p] and p != j]

                if factors_j:
                    # Compute the product of parameters in the same gate
                    product_term_j = np.prod(factors_j)
                    second_term.append((product_term_j, k))  # Append (product, gate index k)
                else:
                    # If no other parameters, append 1
                    second_term.append((1, k))

            # Compute the scalar product between the two terms
            T_ij = 0
            for term_i, n in first_term:
                for term_j, k in second_term:
                    # Multiply the terms and map them to S[n, k] (interaction between gates n and k)
                    T_ij += term_i * term_j * S[n, k]  # S[n, k] refers to the base matrix value

            # Store the result in the T matrix
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
        current_T = np.asmatrix(np.zeros((p + 1, p + 1)))
        current_params = copy.copy(indep_params)
        current_params.append(p)

        for j in range(len(current_params)):
            for k in range(len(current_params)):
                current_T[j, k] = T[current_params[j], current_params[k]]

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



def IDEA(circuit, tol=10 ** (-10), n_points=1000, verbose=True):
    """
    Perform Integrated Dimensional Expressivity Analysis (DEA) using the computed T matrix.

    Parameters:
    - circuit: The quantum circuit being analyzed.
    - tol: Tolerance for considering eigenvalues.
    - n_points: Number of random points to sample.
    - verbose: If True, print progress information.

    Returns:
    - A tuple with:
        - A boolean indicating if the parameters are independent.
        - The set of independent parameters if they exist.
    """
    # Calculate the phi_dict and the phi_values associated to the circuit
    phi_dict, phi_values = compute_phi_dict(circuit)

    # Number of parameters (phi) in the circuit
    n_params = len(phi_dict)
    independent_at_point = []

    if verbose:
        print("Running DEA with T matrix")

    for idx in tqdm.tqdm(range(n_points)):
        # Choose random parameter values for theta (could be phi_values passed)
        theta = 2 * np.pi * np.random.random(n_params)

        # Compute the S matrix for the given circuit and theta values
        S_matrix = compute_S_matrix(circuit, n_params, theta)

        # Compute the T matrix using the computed S matrix, phi_dict, and phi_values
        T = compute_T_matrix(S_matrix, phi_dict, phi_values)

        # Call the modularized function to find independent parameters
        independent_at_point = find_independent_parameters(T, n_params, tol, independent_at_point, theta)

    # Return the result: if only one set of independent parameters was found, return it
    if len(independent_at_point) == 1:
        return True, independent_at_point[0][0]
    return False, independent_at_point



