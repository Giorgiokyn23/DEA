"""
This is an example script running DEA on the circuit

C(\theta) = R_Z(\theta_1) R_X(\theta_0) |0>

"""

from DEA import Circuit, IDEA

C = Circuit(4) 
C.apply(C.H(0))
C.apply(C.H(1))
C.apply(C.H(2))
C.apply(C.H(3))    
C.apply(C.RZ(0)) 
C.apply(C.RZ(1))  
C.apply(C.RZ(2)) 
C.apply(C.RZ(3))
C.apply(C.CNOT(0,1))
C.apply(C.RZ(1))  
C.apply(C.CNOT(0,1))
C.apply(C.CNOT(1,2))
C.apply(C.RZ(2))  
C.apply(C.CNOT(1,2))
C.apply(C.CNOT(2,3))
C.apply(C.RZ(3))  
C.apply(C.CNOT(2,3))
C.apply(C.CNOT(0,3))
C.apply(C.RZ(3))  
C.apply(C.CNOT(0,3))
C.apply(C.CNOT(1,3))
C.apply(C.RZ(3))  
C.apply(C.CNOT(1,3))

# Initialize the phi_dict for this circuit
phi_dict = {
        0: [0, 4, 7],
        1: [1, 4, 5, 8],
        2: [2, 5, 6],
        3: [3, 6, 7, 8]
}


unique1, indep_params1 = IDEA(C,phi_dict,n_points=100,tol=10**(-13))
print(unique1,indep_params1)
if unique1:
    print("The list of independent parameters is "+str(indep_params1)+".")
