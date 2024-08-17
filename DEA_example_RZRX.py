"""
This is an example script running DEA on the circuit

C(\theta) = R_Z(\theta_1) R_X(\theta_0) |0>

"""

import DEA

C = DEA.Circuit(1) # 1 qubit, standard initialization is |0...0>
C.apply(C.RX(0)) # apply X rotation gate to qubit 0
C.apply(C.RZ(0)) # apply Z rotation gate to qubit 0

unique, indep_params = DEA.DEA(C)
if unique:
    print("The list of independent parameters is "+str(indep_params)+".")

unique1, indep_params1 = DEA.DEA_with_T(C)
if unique1:
    print("The list of independent parameters is "+str(indep_params1)+".")
