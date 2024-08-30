# author : pradeep kumar 
# HHL Negativity

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
# import matplotlib.pyplot as plt
# import matplotlib
# from matplotlib import cm
# from mpl_toolkits.mplot3d import axes3d

def PSI_1(b_0,eta_1=1,eta_2=2):
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
    # psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    # psi_a_2 = basis(2,1)
    # psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    # psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    PSI = b_m*tensor(psi_eta_1,psi_c_1) + b_p*tensor(psi_eta_2,psi_c_2)
    # print(PSI)
    return PSI


def PSI_2(b_0,eta_1=1,eta_2=2):
    print(b_0)
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
    for i in range(len_lamda):
        if psi_eta_1.shape[0]==1:
            psi_eta_1 = basis(2,int(eta1_bin[len_lamda-i-1]))
        else : 
            psi_eta_1 = tensor(basis(2,int(eta1_bin[len_lamda-i-1])),psi_eta_1)
        if psi_eta_2.shape[0]==1:
            psi_eta_2 = basis(2,int(eta2_bin[len_lamda-i-1]))
        else :
            psi_eta_2 = tensor(basis(2,int(eta2_bin[len_lamda-i-1])),psi_eta_2)

    # print(psi_eta_1,psi_eta_2)
    # psi_eta_1 = Qobj(psi_eta_1.full().flatten().reshape(psi_eta_1.shape))
    # psi_eta_2 = Qobj(psi_eta_2.full().flatten().reshape(psi_eta_2.shape))
    # print(psi_eta_1,psi_eta_2)

    # return psi_eta_1
    



    # psi_a_1 = basis(2,0)
    # psi_b_1 = basis(2,1)
    # psi_eta_1 = tensor(basis)
    # psi_eta_2 = 
    # print(psi_eta_1)
    # print(psi_eta_2)
    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)
    # print(psi_c_1,psi_d_1)
    # psi_a_2 = basis(2,1)
    # psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)
    # print(psi_c_2,psi_d_2)

    PSI = b_m*tensor(psi_eta_1,psi_c_1,psi_d_1) + b_p*tensor(psi_eta_2,psi_c_2,psi_d_2)
    # print(PSI)
    return PSI

# print(GGM.all_states_matrix(PSI_2(0.3)))

def PSI_3(b_0,eta_1=1,eta_2=2):
    print(b_0)
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    eta1_bin = bin(eta_1)[2:]
    eta2_bin = bin(eta_2)[2:]

    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)
    # print(psi_c_1,psi_d_1)

    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)
    # print(psi_c_2,psi_d_2)

    # alpha = b_m*tensor(psi_c_1,psi_d_1)
    # beta = b_p*tensor(psi_c_2,psi_d_2)
    # # print(alpha, beta)
    # PSI = alpha + beta
    PSI = b_m*tensor(psi_c_1,psi_d_1) + b_p*tensor(psi_c_2,psi_d_2)
    # print(PSI)
    return PSI

def PSI_GGM(b_0,eta_1=1,eta_2=2):
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

    PSI = b_m*tensor(psi_eta_1,psi_c_1,psi_d_1) + b_p*tensor(psi_eta_2,psi_c_2,psi_d_2)
    # print(PSI)
    return PSI

def PSI_2_GGM_disorder(b_0,eta_1=1,eta_2=2,w = 0.1, u = 0.2):
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
    theta1 = np.arcsin(C/eta_1)
    psi_d_1 = np.cos(theta1+w)*basis(2,0) + np.sin(theta1+w)*basis(2,1)

    # psi_a_2 = basis(2,1)
    # psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    theta2 = np.arcsin(C/eta_2)
    psi_d_2 = np.cos(theta2+u)*basis(2,0) + np.sin(theta2+u)*basis(2,1)

    PSI = b_m*tensor(psi_eta_1,psi_c_1,psi_d_1) + b_p*tensor(psi_eta_2,psi_c_2,psi_d_2)
    # print(PSI)
    return PSI

def New_GGM(psi):
    sub_systems = [i for i in range(len(psi.dims[0]))]
    max_lamda_sq = 0
    for i in sub_systems :
        sub_sys = sub_systems[:i] + sub_systems[i+1:]
        # print(sub_sys)
        rho = psi.ptrace(sub_sys)
        eigvals, eigvecs = rho.eigenstates()
        # print(eigvals)
        max_eigs = max(eigvals)
        if max_lamda_sq < max_eigs :
            max_lamda_sq = max_eigs
    ggm = 1 - max_lamda_sq
    return ggm

def max_lambda(psi):
    sub_systems = [i for i in range(len(psi.dims[0]))]
    max_lamda_sq = 0
    for i in sub_systems :
        sub_sys = sub_systems[:i] + sub_systems[i+1:]
        # print(sub_sys)
        rho = psi.ptrace(sub_sys)
        eigvals, eigvecs = rho.eigenstates()
        print(sub_sys,eigvals)
        max_eigs = max(eigvals)
        if max_lamda_sq < max_eigs :
            max_lamda_sq = max_eigs
    return max_eigs
    # ggm = 1 - max_lamda_sq
    # return ggm



def PSI_2_2(b_0) : 

    # b_0 = 0.2
    print(b_0)
    b_1 = np.sqrt(1-b_0**2)
    print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    eta_1 = 1.0
    eta_2 = 2.0
    C = 0.736

    psi_a_1 = basis(2,0)
    psi_b_1 = basis(2,1)
    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_a_2 = basis(2,1)
    psi_b_2 = basis(2,0)
    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    PSI = b_m*tensor(psi_a_1,psi_b_1,psi_c_1,psi_d_1) + b_p*tensor(psi_a_2,psi_b_2,psi_c_2,psi_d_2)
    # print(PSI)
    return PSI

# rho_ab = PSI.ptrace([0,2])
# rho_ac = PSI.ptrace([0,3])
# print(rho_ab)
# print(rho_ac)


def get_negativity_b0(index_list):
    N_abc = []
    b_0_square =[]
    b_0 = [i/1000.0 for i in range(0,1001) ]

    for j in range(len(b_0)):
        PSI = PSI_2(b_0[j])
        b_0_square.append(b_0[j]**2)
        # N_abc.append(qutip.negativity(PSI.ptrace([0,1,2]),0))
        N_abc.append(qutip.negativity(PSI.ptrace(index_list),0))
    
    return (b_0_square,N_abc)

def plot_negativity(index_list):
    b_0_square,N_abc = get_negativity_b0(index_list)

    plt.plot(b_0_square,N_abc)
    plt.show()

def custom_negativity(PSI, sys , sub_sys ):
    if len(sys) == int(np.log2(PSI.shape[0])):
        rho = PSI*PSI.dag()
    else :
        rho = PSI.ptrace(sys)
    rho_T = qutip.partial_transpose(rho,sub_sys)
    # sys should be the index which remains after partial trace
    # sub_sys is the list which indicate which should be transposed 
    # ex. [1,0,1] this means first and last index should be transposed
    def get_neg_eigh_sum(eigens):
        s = 0
        for i in eigens:
            if i < 0 :
                # print(i)
                s = s+i
        return -1*s
    return get_neg_eigh_sum(LA.eigh(rho_T)[0])



def negativity_system(PSI):
    qubits = int(np.log2(PSI.shape[0]))
    sys = [i for i in range(qubits)]
    subsys = [0 for i in range(qubits)]
    subsys[-1] = 1
    return custom_negativity(PSI, sys, subsys)

def get_negativity(eta_1=1,eta_2=2):
    b_0 = [i/1000.0 for i in range(1001)]
    negativity = []
    b_0_sq = [i**2 for i in b_0]
    for i in range(len(b_0)):
        PSI = PSI_2(b_0[i],eta_1, eta_2)
        qubits = int(np.log2(PSI.shape[0]))
        sub_sys = [0 for i in range(qubits)]
        sub_sys[-1] = 1
        negativity.append(custom_negativity(PSI,[i for i in range(qubits)],sub_sys))
    return (b_0_sq,negativity)
    

def print_all_negativities(PSI):
    qubits = int(np.log2(PSI.shape[0]))
    qubit_index = [i for i in range(qubits)]
    # print(qubits)
    for k in range(2,qubits+1):
        comb = list(combinations(qubit_index,k))
        # print(comb) 
        for i in comb: 
            sub_sys = list(product(range(2),repeat=len(i)))
            for j in sub_sys:
                # sub_sys = product
                print(i, j, custom_negativity(PSI, i , j))
                # pass
    return


def plot_full_neg():
    N = []
    b_0_square =[]
    b_0 = [i/1000.0 for i in range(0,1001)]
    def get_neg_eigh_sum(eigens):
        s = 0
        for i in eigens:
            if i < 0 :
                s = s+i
        return -1*s
    for i in b_0:
        PSI = PSI_2(i)
        rho_psi = PSI*PSI.dag()
        mask_len = int(np.log2(PSI.shape[0]))
        mask = [0 for i in range(mask_len)]
        mask[0] = 1
        rho_pt = qutip.partial_transpose(rho_psi,mask)
        neg = get_neg_eigh_sum(LA.eigh(rho_pt)[0])
        N.append(neg)
        b_0_square.append(i**2)
    
    plt.plot(b_0_square,N)
    plt.show()



def GGM(PSI):


    lamda_PSI = []
    shape = PSI.shape[0]
    if np.log2(shape)%2 == 0 :
        ul = int(np.log2(shape)/2.0)
    else :
        ul = int(np.log2(shape)/2.0) + 1
    for i in range(1,ul+1):
        cut_iterations = list(combinations([l for l in range(int(np.log2(shape)))],i))
        # print(cut_iterations)
        for cut in cut_iterations : 
            # print(cut)
            PSI_new = GGM_module.all_states_matrix(PSI, cut)
            eigens = LA.svd(PSI_new)[1]
            eigens_sq = [i**2 for i in eigens]
            # print(max(eigens_sq))
            lamda_PSI.append(max(eigens_sq))

    max_GGM = max(lamda_PSI)

    ggm = 1 - max_GGM
    return ggm


def get_GGM(eta_1=1, eta_2 =2):
    b_0 = [i/1000.0 for i in range(1001)]
    ggm = []
    b_0_sq = [i**2 for i in b_0]
    for i in range(len(b_0)):
        PSI = PSI_2(b_0[i],eta_1, eta_2)
        ggm.append(GGM(PSI))
    return (b_0_sq,ggm)

def plot_GGM():
    b_0_sq, ggm = get_GGM()
    plt.plot(b_0_sq, ggm)
    plt.show()


def plot_both(index_list):
    #plots GGM and negativity
    ggm = get_GGM()[1]
    b_0,neg = get_negativity_b0(index_list)
    plt.plot(b_0,neg)
    plt.plot(b_0,ggm)
    plt.show()
