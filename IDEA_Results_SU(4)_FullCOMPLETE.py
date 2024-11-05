from DEA import Circuit, DEA, IDEA


C = Circuit(4)
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

#Beginning of the last layer
C.apply(C.RZ(2))       # RZ gate on qubit1
C.apply(C.CNOT(1,2))  # CNOT between qubit0 and qubit1
C.apply(C.RZ(1))       # RZ gate on qubit0
C.apply(C.RX(1))       # RX gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit1
C.apply(C.RX(2))       # RX gate on qubit1
C.apply(C.RZ(2))       # RZ gate on qubit1
C.apply(C.CNOT(1, 2))  # CNOT between qubit0 and qubit1
C.apply(C.RX(1))       # RX gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit1
C.apply(C.CNOT(1, 2))  # CNOT between qubit0 and qubit1
C.apply(C.RZ(1))       # RZ gate on qubit0
C.apply(C.RX(1))       # RX gate on qubit0
C.apply(C.RZ(1))       # RZ gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit1
C.apply(C.RX(2))       # RX gate on qubit1
C.apply(C.RZ(2))
#Adding from 15 parameters
C.apply(C.RZ(0))       # RZ gate on qubit0
C.apply(C.RY(0))       # RX gate on qubit0
C.apply(C.RZ(0))
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RY(2))       # RX gate on qubit0
C.apply(C.RZ(2))
C.apply(C.CNOT(0, 2))
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RY(2))       # RX gate on qubit0
C.apply(C.RZ(2))
C.apply(C.RZ(1))       # RZ gate on qubit0
C.apply(C.RY(1))       # RX gate on qubit0
C.apply(C.RZ(1))
C.apply(C.RZ(3))       # RZ gate on qubit0
C.apply(C.RY(3))       # RX gate on qubit0
C.apply(C.RZ(3))
C.apply(C.CNOT(1, 3))
C.apply(C.RZ(3))       # RZ gate on qubit0
C.apply(C.RY(3))       # RX gate on qubit0
C.apply(C.RZ(3))

# Last Single Convolutional Layer

C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.CNOT(2,3))  # CNOT between qubit0 and qubit1
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RX(2))       # RX gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.RX(3))       # RX gate on qubit1
C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.CNOT(2, 3))  # CNOT between qubit0 and qubit1
C.apply(C.RX(2))       # RX gate on qubit0
C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.CNOT(2, 3))  # CNOT between qubit0 and qubit1
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RX(2))       # RX gate on qubit0
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RZ(3))       # RZ gate on qubit1
C.apply(C.RX(3))       # RX gate on qubit1
C.apply(C.RZ(3))

#Last LAyer
C.apply(C.RZ(2))       # RZ gate on qubit0
C.apply(C.RY(2))       # RX gate on qubit0
C.apply(C.RZ(2))
C.apply(C.RZ(3))       # RZ gate on qubit0
C.apply(C.RY(3))       # RX gate on qubit0
C.apply(C.RZ(3))
C.apply(C.CNOT(2, 3))
C.apply(C.RZ(3))       # RZ gate on qubit0
C.apply(C.RY(3))       # RX gate on qubit0
C.apply(C.RZ(3))


phi_dict = {
    0: [0,1,30],
    1: [2,8,31],
    2: [3,9,32],
    3: [4,10,33],
    4: [5,11,34],
    5: [6,12,35],
    6: [7,13,36],
    7: [14,16,37],
    8: [15,17,38],
    9: [18,24,39],
    10: [19,25,40],
    11: [20,26,41],
    12: [21,27,42],
    13: [22,28,43],
    14: [23,29,44],
    15: [45, 54],
    16: [46, 55],
    17: [47, 56],
    18: [48, 57],
    19: [49, 58],
    20: [50, 59],
    21: [51, 60],
    22: [52, 61],
    23: [53, 62],
    24: [63],
    25: [64],
    26: [65],
    27: [66],
    28: [67],
    29: [68],
    30: [69],
    31: [70],
    32: [71],
    33: [72],
    34: [73],
    35: [74],
    36: [75],
    37: [76],
    38: [77],
    39: [78],
    40: [79],
    41: [80],
    42: [81],
    43: [82],
    44: [83],
    45: [84],
    46: [85],
    47: [86],
}

unique1, indep_params1 = IDEA(C,phi_dict,n_points=2,tol=10**(-13))
print(unique1,indep_params1)
if unique1:
    print("The list of independent parameters is "+str(indep_params1)+".")