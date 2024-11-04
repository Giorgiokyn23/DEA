from DEA import Circuit, DEA, IDEA


C = Circuit(2)
C.reset_random()
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RZ(3))
C.apply(C.CNOT(0, 1))  # CNOT between qubit0 and qubit1
C.apply(C.CNOT(2,3))
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RX(0))       # RX gate on qubit0
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RX(1))       # RX gate on qubit1
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RX(2))       # RX gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.RX(3))       # RX gate on qubit1
C.apply(C.RZ(3))
C.apply(C.CNOT(0, 1))
C.apply(C.CNOT(2, 3))

C.apply(C.RX(0))       # RX gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RX(2))       # RX gate on qubit0
C.apply(C.RZ(3))
C.apply(C.CNOT(0, 1))
C.apply(C.CNOT(2, 3))

C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RX(0))       # RX gate on qubit0
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit1
C.apply(C.RX(1))       # RX gate on qubit1
C.apply(C.RZ(1))

C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RX(2))       # RX gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.RX(3))       # RX gate on qubit1
C.apply(C.RZ(3))

phi_dict = {
    0: [0,1],
    1: [2,8],
    2: [3,9],
    3: [4,10],
    4: [5,11],
    5: [6,12],
    6: [7,13],
    7: [14,16],
    8: [15,17],
    9: [18,24],
    10: [19,25],
    11: [20,26],
    12: [21,27],
    13: [22,28],
    14: [23,29],
}

unique1, indep_params1 = IDEA(C,phi_dict,n_points=100,tol=10**(-13))
print(unique1,indep_params1)
if unique1:
    print("The list of independent parameters is "+str(indep_params1)+".")