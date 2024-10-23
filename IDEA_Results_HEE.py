from DEA import Circuit, IDEA

C = Circuit(4)
C.apply(C.RX(0))
C.apply(C.RX(1))
C.apply(C.RX(2))
C.apply(C.RX(3))
C.apply(C.CNOT(0,1))
C.apply(C.CNOT(1,2))
C.apply(C.CNOT(2,3))
C.apply(C.RX(0))
C.apply(C.RX(1))
C.apply(C.RX(2))
C.apply(C.RX(3))
C.apply(C.CNOT(0,1))
C.apply(C.CNOT(1,2))
C.apply(C.CNOT(2,3))
C.apply(C.RX(0))
C.apply(C.RX(1))
C.apply(C.RX(2))
C.apply(C.RX(3))
C.apply(C.CNOT(0,1))
C.apply(C.CNOT(1,2))
C.apply(C.CNOT(2,3))
# running with 4 parameters there are 4 independent parameters

# Initialize the phi_dict for this circuit
phi_dict = {
        0: [0,4,8],
        1: [1,5,9],
        2: [2,6,10],
        3: [3,7,11]
}


unique1, indep_params1 = IDEA(C,phi_dict,n_points=100,tol=10**(-13))
print(unique1,indep_params1)
if unique1:
    print("The list of independent parameters is "+str(indep_params1)+".")