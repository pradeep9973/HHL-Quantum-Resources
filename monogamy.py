#author: pradeep kumar
# plotting monogamy


import GGM_module
from itertools import combinations
from qutip.states import *
from qutip.tensor import *
import matplotlib as mpl
from qutip.qobj import Qobj, isket
import qutip
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import HHL_negativity as hhl

psi_list = [hhl.PSI_2(0.2,i,100) for i in range(2,100)]
def Monogamy(psi):
    qubits = int(np.log2(psi.shape[0]))
    sys = [i for i in range(qubits)]
    sub_sys = [0 for i in range(qubits)]
    sub_sys[-1] = 1
    cn_1 = hhl.custom_negativity(psi,sys,sub_sys)
    cn_2 = hhl.custom_negativity(psi,sys[:-1],sub_sys[1:])
    mono = cn_2**2 - cn_1**2
    return (cn_1,cn_2,mono)

ratio = ratio = [i / 100.0 for i in range(2, 100)]

mono = []
cn_1 = []
cn_2 = []
counter = 0
for i in psi_list:
    # mono.append(Monogamy(i))
    a,b,c = Monogamy(i)
    mono.append(c)
    cn_1.append(a**2)
    cn_2.append(b**2)
    print(counter)
    counter += 1

plt.plot(ratio,mono,label = "mono")
plt.plot(ratio,cn_1,label = "cn_1")
plt.plot(ratio, cn_2,label = "cn_2")
plt.legend()

plt.show()

