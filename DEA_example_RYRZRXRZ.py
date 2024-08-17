"""
This is an example script running DEA on the circuit

C(\theta) = R_Y(\theta_3) R_Z(\theta_2) R_X(\theta_1) R_Z(\theta_0) |0>

"""

import DEA

C = DEA.Circuit(1) # 1 qubit, standard initialization is |0...0>
C.apply(C.RZ(0)) # apply Z rotation gate to qubit 0
C.apply(C.RX(0)) # apply X rotation gate to qubit 0
C.apply(C.RZ(0)) # apply Z rotation gate to qubit 0
C.apply(C.RY(0)) # apply Y rotation gate to qubit 0

unique, indep_params = DEA.DEA(C)
if unique:
    print("The list of independent parameters is "+str(indep_params)+".")
else:
    print(indep_params)
