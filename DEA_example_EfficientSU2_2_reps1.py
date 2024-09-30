"""
This is an example script running DEA on the QISKIT circuit EfficientSU2(2,reps=1) 

|0> - RY(\theta_0) - RZ(\theta_2) - CNOT control - RY(\theta_4) - RZ(\theta_6)
|0> - RY(\theta_1) - RZ(\theta_3) - CNOT target  - RY(\theta_5) - RZ(\theta_7)

"""
import numpy as np
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

phi_dict,phi_values=DEA.compute_phi_dict(C)
print(phi_dict, phi_values)

#S= DEA.compute_S_matrix(C,8,1.234)
#print(S)

S = np.array([
    [0.9, 0.1, 0.2, 0.0, 0.3, 0.0, 0.4, 0.0],
    [0.1, 0.8, 0.3, 0.1, 0.0, 0.2, 0.0, 0.4],
    [0.2, 0.3, 0.7, 0.0, 0.4, 0.0, 0.1, 0.0],
    [0.0, 0.1, 0.0, 0.6, 0.0, 0.3, 0.0, 0.2],
    [0.3, 0.0, 0.4, 0.0, 0.5, 0.0, 0.2, 0.0],
    [0.0, 0.2, 0.0, 0.3, 0.0, 0.4, 0.0, 0.1],
    [0.4, 0.0, 0.1, 0.0, 0.2, 0.0, 0.6, 0.0],
    [0.0, 0.4, 0.0, 0.2, 0.0, 0.1, 0.0, 0.7],
])

# Step 4: Compute the T matrix based on S, phi_dict, and phi_values
T = DEA.compute_T_matrix(S, phi_dict, phi_values)

# Step 5: Print the T matrix result
print("T matrix calculated:")
print(T)
#print(C._param_to_gate)



unique, indep_params = DEA.IDEA(C,n_points=100)
if unique:
    print("The list of independent parameters with S is "+str(indep_params)+".")
else:
    print(indep_params)
"""
C above is a base circuit >> S matrix in DEA
base parameters >> theta

now, we have new parameters >> phi
phi parameters get mapped to theta parameters
phi_0 >> theta_1 and theta_3
phi_1 >> theta_0 and theta_2

phi_dict = {0:[1,3], 1:[0,2]}

T[0,1] = S[1,0]+S[1,2]+S[3,0]+S[3,2]

---------------

Let's say we make EfficientSU2 translationally invariant >> gates at the same time have the same parameter

|0> - RY(theta0) - RZ(theta2) - C - RY(theta4) - RZ(theta6)
|0> - RY(theta1) - RZ(theta3) - X - RY(theta5) - RZ(theta7)

same time = same parameter: new circuit

|0> - RY(phi0) - RZ(phi1) - C - RY(phi2) - RZ(phi3)
|0> - RY(phi0) - RZ(phi1) - X - RY(phi2) - RZ(phi3)

phi = {0:[0,1], 1:[2,3], 2:[4,5], 3:[6,7]}

-----------
expected input : base circuit + a data structure that tells me how trainable parameters are mapped to base circuit parameters (phi_dict)

phi_dict = {0:[0,2], 1:[1,2]}

|0> - RY(phi0) - RZ(phi1) - RX(phi0*phi1)

T[0,1] = (d_0+phi1*d_2),(d_1+phi0*d_2) = S[0,1] + phi1*S[2,1] + phi0*S[0,2] + phi0*phi1*S[2,2]

-----

parser possibiity:

rather than phi: trainable parameters to base parameters

use 2 objects
object 1: list base parameters that contain trainable parameters, 
e.g. for 
|0> - RY(phi0) - RZ(phi1) - RX(phi0*phi)
O1 = {0:[0],1:[1],2:[0,1]}

object 2: expressions of trainable parameters
O2 = {0:"p0", 1:"p1", 2:"p0*p1"}

3:"2*p0+3*p1"
|0> - RY(phi0) - RZ(phi1) - RX(phi0*phi) - RY(2*phi0+3*phi1)

"""



"""
unique, indep_params = DEA.DEA(C,n_points=100)
if unique:
    print("The list of independent parameters with S is "+str(indep_params)+".")
else:
    print(indep_params)

unique1, indep_params1 = DEA.DEA_with_T(C,n_points=100)
if unique1:
    print("The list of independent parameters with T is "+str(indep_params1)+".")
else:
    print(indep_params1)

"""
