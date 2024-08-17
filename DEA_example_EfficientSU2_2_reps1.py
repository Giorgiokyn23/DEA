"""
This is an example script running DEA on the QISKIT circuit EfficientSU2(2,reps=1) 

|0> - RY(\theta_0) - RZ(\theta_2) - CNOT control - RY(\theta_4) - RZ(\theta_6)
|0> - RY(\theta_1) - RZ(\theta_3) - CNOT target  - RY(\theta_5) - RZ(\theta_7)

"""

import DEA

C = DEA.Circuit(2) # 2 qubit, standard initialization is |0...0>
C.apply(C.RY(0))         # apply Y rotation gate to qubit 0
C.apply(C.RY(1))         # apply Y rotation gate to qubit 1
C.apply(C.RZ(0))         # apply Z rotation gate to qubit 0
C.apply(C.RZ(1))         # apply Z rotation gate to qubit 1
C.apply(C.CNOT(0,1)) # apply CNOT with control 0 and target 1
C.apply(C.RY(0))         # apply Y rotation gate to qubit 0
C.apply(C.RY(1))         # apply Y rotation gate to qubit 1
C.apply(C.RZ(0))         # apply Z rotation gate to qubit 0
C.apply(C.RZ(1))         # apply Z rotation gate to qubit 1

unique, indep_params = DEA.DEA(C,n_points=100)
if unique:
    print("The list of independent parameters is "+str(indep_params)+".")
else:
    print(indep_params)
