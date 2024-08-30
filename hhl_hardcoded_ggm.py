#hhl_hardcoded_ggm.py


# import GGM_module
# from itertools import combinations
# from itertools import product
from qutip.states import *
from qutip.tensor import *
# import matplotlib as mpl
from qutip.qobj import Qobj, isket
import qutip
import numpy as np
from numpy import linalg as LA
import random
import HHL_negativity
import matplotlib.pyplot as plt
import time 
import hhl_entanglement





def PSI_2_GGM_disorder(b_0,w = 0.1, u = 0.2):
    # print(b_0)
    
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    eta_1 = 1.0
    eta_2 = 2.0
    C = 0.736

    psi_eta_1 = qutip.basis(4,1)  #tensor(basis(2,0),basis(2,1))
    psi_eta_2 = qutip.basis(4,2)  #tensor(basis(2,1),basis(2,0))
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

def HC_GGM(psi):
    # sub_systems = [i for i in range(len(psi.dims[0]))]
    sub_systems = [0,1,2]
    max_lamda_sq = 0
    for i in sub_systems :
        rho = psi.ptrace([i])
        eigvals, eigvecs = LA.eig(rho) #rho.eigenstates()
        # print(eigvals)
        # eigens = eigens.join(eigvals)
        max_eigs = max(eigvals)
        if max_lamda_sq < max_eigs :
            max_lamda_sq = max_eigs
    ggm = 1 - max_lamda_sq
    return ggm


def ggm_disorder(b_0,sigma = 0.2):
    rn1 = random.gauss(0,sigma)
    rn2 = random.gauss(0,sigma)
    psi = PSI_2_GGM_disorder(b_0,rn1,rn2)
    return HC_GGM(psi), hhl_entanglement.get_coherence(psi.ptrace([2])), hhl_entanglement.get_lneg(psi,[0,1],[0,1])

def get_ggm_disorder_data(sigma=0.2):
    
    number_of_runs = 100

    number_of_steps = 100

    ggm = []
    C_r =[]
    b0_sq = [i/(1.0*number_of_steps) for i in range(number_of_steps+1)]

    for i in b0_sq:
        sum_ggm = 0
        sum_coherence = 0
        for j in range(number_of_runs):
            a,b,c = ggm_disorder(np.sqrt(i),sigma)
            sum_ggm += a
            sum_coherence += b
        ggm.append(sum_ggm/number_of_runs)
        C_r.append(sum_coherence/number_of_runs)
        print(i,"*****************")
    return b0_sq,ggm,C_r


def avg_entanglement(b_0,sigma=0.2):
    start = time.time()
    number_of_runs = 10000

    # number_of_steps = 100

    # ggm = []
    # C_r =[]
    # b0_sq = #[i/(1.0*number_of_steps) for i in range(number_of_steps+1)]


    sum_ggm = 0
    sum_coherence = 0
    sum_LN = 0
    for j in range(number_of_runs):
        a,b,c = ggm_disorder(b_0,sigma)
        sum_ggm += a
        sum_coherence += b
        sum_LN += c
    avg_ggm = (sum_ggm/number_of_runs)
    avg_C_r = (sum_coherence/number_of_runs)
    avg_LN = (sum_LN/number_of_runs)
    # print(i,"*****************")
    print(b_0,b_0**2, avg_ggm,avg_C_r,avg_LN)
    # return ggm,C_r
    print(time.time()-start)

def plot_disorder():
    start = time.time()
    b0_sq,ggm = get_ggm_disorder_data()
    plt.plot(b0_sq,ggm)
    end = time.time()
    print(end-start)
    plt.show()
    
    
def get_error(psi3_og, psi3_d):

    # return abs(np.dot(a1.conj(),a2)/(abs(a1)*abs(a2)))
    # return abs(np.dot(a1.conj(),a2) - np.dot(a2.conj(),a2))
    a1 = np.array([psi3_og[1][0][0],psi3_og[3][0][0]])
    a2 = np.array([psi3_d[1][0][0],psi3_d[3][0][0]])

    # print(np.linalg.norm(a1),np.linalg.norm(a2))
    return np.linalg.norm(a1-a2)/np.linalg.norm(a1)

def PSI_3_HC(b_0):
    
    eta_1=1
    eta_2=2  

    # print("b_0 : ", b_0)
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    # eta1_bin = bin(eta_1)[2:]
    # eta2_bin = bin(eta_2)[2:]

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

def PSI_3_disorder_HC(b_0,w=0.1,u=0.2):
    # print(b_0)
    eta_1=1
    eta_2=2
    b_1 = np.sqrt(1-b_0**2)
    # print(b_1)
    b_m = (b_0-b_1)/np.sqrt(2)
    b_p = (b_0+b_1)/np.sqrt(2)


    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    # eta1_bin = bin(eta_1)[2:]
    # eta2_bin = bin(eta_2)[2:]

    psi_c_1 = (basis(2,0) - basis(2,1)).unit()
    theta1 = np.arcsin(C/eta_1)
    # print(theta1, np.sin(theta1), np.sin(theta1 + w))
    # psi_d_1_og = np.cos(theta1)*basis(2,0) + np.sin(theta1)*basis(2,1)
    psi_d_1 = np.cos(theta1+w)*basis(2,0) + np.sin(theta1+w)*basis(2,1)
    # print(psi_d_1_og, psi_d_1)

    psi_c_2 = (basis(2,0) + basis(2,1)).unit()
    theta2 = np.arcsin(C/eta_2)
    psi_d_2 = np.cos(theta2+u)*basis(2,0) + np.sin(theta2+u)*basis(2,1)


    PSI = b_m*tensor(psi_c_1,psi_d_1) + b_p*tensor(psi_c_2,psi_d_2)
    # print(PSI)
    return PSI


def get_avg_error(b_0,sigma):
    start = time.time()
    runs = 10000

    sum_error = 0
    for i in range(runs):
        rn1 = random.gauss(0,sigma)
        rn2 = random.gauss(0,sigma)
        sum_error += get_error(PSI_3_HC(b_0),PSI_3_disorder_HC(b_0,rn1,rn2))
    print(time.time() - start)
    return sum_error/runs

        



