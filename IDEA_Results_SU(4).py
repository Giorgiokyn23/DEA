from DEA import Circuit, DEA, IDEA


C = Circuit(2)
C.reset_random()
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.CNOT(0, 1))  # CNOT between qubit0 and qubit1
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RX(0))       # RX gate on qubit0
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RX(1))       # RX gate on qubit1
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.CNOT(0, 1))  # CNOT between qubit0 and qubit1
C.apply(C.RX(0))       # RX gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.CNOT(0, 1))  # CNOT between qubit0 and qubit1
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RX(0))       # RX gate on qubit0
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RX(1))       # RX gate on qubit1
C.apply(C.RZ(1))

phi_dict = {
    0: [0],
    1: [1],
    2: [2],
    3: [3],
    4: [4],
    5: [5],
    6: [6],
    7: [7],
    8: [8],
    9: [9],
    10: [10],
    11: [11],
    12: [12],
    13: [13],
    14: [14],
}

unique1, indep_params1 = IDEA(C,phi_dict,n_points=100,tol=10**(-13))
print(unique1,indep_params1)
if unique1:
    print("The list of independent parameters is "+str(indep_params1)+".")
