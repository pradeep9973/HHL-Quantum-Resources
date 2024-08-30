#ggm for different partitions

from math import ceil
import GGM_module
from itertools import combinations
from itertools import product
from qutip.states import *
from qutip.tensor import *
# import matplotlib as mpl
from qutip.qobj import Qobj, isket
import qutip
import numpy as np
from numpy import linalg as LA
import HHL_negativity



## this way of calculating GGM/coherence is wrong
def PSI_GGM_p1(b_0,eta_1=1,eta_2=2):
    # partition lambda and U
    # print(b_0) 

    
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    eta1_bin = bin(eta_1)[2:]
    eta2_bin = bin(eta_2)[2:]
    # print(eta1_bin,eta2_bin)
    psi_eta_1 = tensor(basis(1))
    psi_eta_2 = tensor(basis(1))
    len_lamda = max(len(eta1_bin),len(eta2_bin))
    if len(eta1_bin) < len(eta2_bin):
        eta1_bin = "0"*(len(eta2_bin) - len(eta1_bin)) + eta1_bin
    elif len(eta2_bin) < len(eta1_bin):
        eta2_bin = "0"*(len(eta1_bin) - len(eta2_bin)) + eta2_bin 
    # print(eta1_bin,eta2_bin)  
    # print(eta1_bin,eta2_bin)  

    len_lamda = len(eta2_bin)

    for i in range(len_lamda):

        psi_eta_1 = tensor(basis(2,int(eta1_bin[len_lamda-i-1])),psi_eta_1)
        psi_eta_2 = tensor(basis(2,int(eta2_bin[len_lamda-i-1])),psi_eta_2)

    # print(psi_eta_1,psi_eta_2)
    psi_eta_1 = Qobj(psi_eta_1.full().flatten().reshape(psi_eta_1.shape))
    psi_eta_2 = Qobj(psi_eta_2.full().flatten().reshape(psi_eta_2.shape))
    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    # psi_a_2 = basis(2,1)
    # psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    PSI = b_m*tensor(psi_eta_1,psi_c_1) + b_p*tensor(psi_eta_2,psi_c_2)
    # print(PSI)
    return PSI


def PSI_GGM_p2(b_0,eta_1=1,eta_2=2):
    # partition lambda and R
    # print(b_0) 

    
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    eta1_bin = bin(eta_1)[2:]
    eta2_bin = bin(eta_2)[2:]
    # print(eta1_bin,eta2_bin)
    psi_eta_1 = tensor(basis(1))
    psi_eta_2 = tensor(basis(1))
    len_lamda = max(len(eta1_bin),len(eta2_bin))
    if len(eta1_bin) < len(eta2_bin):
        eta1_bin = "0"*(len(eta2_bin) - len(eta1_bin)) + eta1_bin
    elif len(eta2_bin) < len(eta1_bin):
        eta2_bin = "0"*(len(eta1_bin) - len(eta2_bin)) + eta2_bin 
    # print(eta1_bin,eta2_bin)  
    # print(eta1_bin,eta2_bin)  

    len_lamda = len(eta2_bin)

    for i in range(len_lamda):

        psi_eta_1 = tensor(basis(2,int(eta1_bin[len_lamda-i-1])),psi_eta_1)
        psi_eta_2 = tensor(basis(2,int(eta2_bin[len_lamda-i-1])),psi_eta_2)

    # print(psi_eta_1,psi_eta_2)
    psi_eta_1 = Qobj(psi_eta_1.full().flatten().reshape(psi_eta_1.shape))
    psi_eta_2 = Qobj(psi_eta_2.full().flatten().reshape(psi_eta_2.shape))
    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    # psi_a_2 = basis(2,1)
    # psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    PSI = b_m*tensor(psi_eta_1,psi_d_1) + b_p*tensor(psi_eta_2,psi_d_2)
    # print(PSI)
    return PSI


def PSI_GGM_p3(b_0,eta_1=1,eta_2=2):
    # partition U and R
    # print(b_0) 

    
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    eta1_bin = bin(eta_1)[2:]
    eta2_bin = bin(eta_2)[2:]
    # print(eta1_bin,eta2_bin)
    psi_eta_1 = tensor(basis(1))
    psi_eta_2 = tensor(basis(1))
    len_lamda = max(len(eta1_bin),len(eta2_bin))
    if len(eta1_bin) < len(eta2_bin):
        eta1_bin = "0"*(len(eta2_bin) - len(eta1_bin)) + eta1_bin
    elif len(eta2_bin) < len(eta1_bin):
        eta2_bin = "0"*(len(eta1_bin) - len(eta2_bin)) + eta2_bin 
    # print(eta1_bin,eta2_bin)  
    # print(eta1_bin,eta2_bin)  

    len_lamda = len(eta2_bin)

    for i in range(len_lamda):

        psi_eta_1 = tensor(basis(2,int(eta1_bin[len_lamda-i-1])),psi_eta_1)
        psi_eta_2 = tensor(basis(2,int(eta2_bin[len_lamda-i-1])),psi_eta_2)

    # print(psi_eta_1,psi_eta_2)
    psi_eta_1 = Qobj(psi_eta_1.full().flatten().reshape(psi_eta_1.shape))
    psi_eta_2 = Qobj(psi_eta_2.full().flatten().reshape(psi_eta_2.shape))
    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    # psi_a_2 = basis(2,1)
    # psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    PSI = b_m*tensor(psi_c_1,psi_d_1) + b_p*tensor(psi_c_2,psi_d_2)
    # print(PSI)
    return PSI


