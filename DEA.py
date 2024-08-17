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


def compute_T_matrix(S,param_gate_map):
    """
   Compute the T matrix by applying the product and chain rule.
   param_gate_map: A dictionary mapping each parameter to the gates it affects.
            For example,{0: [1,5], 1:[2,3,7]} indicates that parameter 0 is in gates 1 and 5 and parameter 1 is in gate 2,3,7

    """

    n_params = len(param_gate_map)
    T = np.asmatrix((np.zeros((n_params,n_params)))

    # Apply product rule accross gates
    for j in range(n_params):
        gates_j = param_gate_map[j] # gates affected by parameter j
        for k in range(n_params):
            gates_k= param_gate_map[k] # gates affected by parameter k
            # Sum the S matrix elements corresponding to gate combinations for j and k