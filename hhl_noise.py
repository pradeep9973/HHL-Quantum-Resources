# hhl with noise

from math import ceil
import HHL_negativity
# import GGM_module
from itertools import combinations
from itertools import product
from qutip.states import *
from qutip.tensor import *
from qutip.qobj import Qobj, isket
import qutip
import numpy as np
from numpy import linalg as LA

def PSI_3d_GGM(b_0=0.5,b_1=0.4,eta_1=1,eta_2=2,eta_3=3):
    # here the eigenvalue is represented by a single vector

    
    b_2 = np.sqrt(1 - b_0**2 - b_1**2)
    # print(b_2)


    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    beta_1 = (b_0 + b_1 + b_2)/np.sqrt(3)
    beta_2 = (b_1 - b_2)/np.sqrt(2)
    beta_3 = (-2*b_0 + b_1 + b_2)/np.sqrt(6)

    

    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    eta1_bin = bin(eta_1)[2:]
    eta2_bin = bin(eta_2)[2:]
    eta3_bin = bin(eta_3)[2:]

    psi_eta_1 = tensor(basis(1))
    psi_eta_2 = tensor(basis(1))
    psi_eta_3 = tensor(basis(1))
    len_lambda = max(len(eta1_bin),len(eta2_bin),len(eta3_bin))

 
    if len(eta1_bin) < len_lambda:
        eta1_bin = "0"*(len_lambda - len(eta1_bin)) + eta1_bin
    if len(eta2_bin) < len_lambda:
        eta2_bin = "0"*(len_lambda - len(eta2_bin)) + eta2_bin
    if len(eta3_bin) < len_lambda:
        eta3_bin = "0"*(len_lambda - len(eta3_bin)) + eta3_bin
    # print(eta1_bin,eta2_bin)  
    # print(eta1_bin,eta2_bin,eta3_bin)  
 
    # len_lamda = len(eta2_bin)

    for i in range(len_lambda):

        psi_eta_1 = tensor(basis(2,int(eta1_bin[len_lambda-i-1])),psi_eta_1)
        psi_eta_2 = tensor(basis(2,int(eta2_bin[len_lambda-i-1])),psi_eta_2)
        psi_eta_3 = tensor(basis(2,int(eta3_bin[len_lambda-i-1])),psi_eta_3)

    # print(psi_eta_1,psi_eta_2)
    psi_eta_1 = Qobj(psi_eta_1.full().flatten().reshape(psi_eta_1.shape))
    psi_eta_2 = Qobj(psi_eta_2.full().flatten().reshape(psi_eta_2.shape))
    psi_eta_3 = Qobj(psi_eta_3.full().flatten().reshape(psi_eta_3.shape))

    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(3,0) + basis(3,1) + basis(3,2) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(3,1) - basis(3,2)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (-2*basis(3,0) + basis(3,1) + basis(3,2)).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    PSI = beta_1*tensor(psi_eta_1,psi_c_1,psi_d_1) + beta_2*tensor(psi_eta_2,psi_c_2,psi_d_2) + beta_3*tensor(psi_eta_3,psi_c_3,psi_d_3)
    # print(PSI)
    return PSI


def PSI_3_3d_GGM(b_0=0.5,b_1=0.4,eta_1=1,eta_2=2,eta_3=3):
    # here the eigenvalue is represented by a single vector

    
    b_2 = np.sqrt(1 - b_0**2 - b_1**2)
    # print(b_2)


    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    beta_1 = (b_0 + b_1 + b_2)/np.sqrt(3)
    beta_2 = (b_1 - b_2)/np.sqrt(2)
    beta_3 = (-2*b_0 + b_1 + b_2)/np.sqrt(6)

    

    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736

    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(3,0) + basis(3,1) + basis(3,2) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(3,1) - basis(3,2)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (-2*basis(3,0) + basis(3,1) + basis(3,2)).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    PSI = beta_1*tensor(psi_c_1,psi_d_1) + beta_2*tensor(psi_c_2,psi_d_2) + beta_3*tensor(psi_c_3,psi_d_3)
    # print(PSI)
    return PSI




def GGM3d(b_0 = 0.4, b_1 = 0.3):
    a = HHL_negativity.New_GGM(PSI_3d_GGM(b_0,b_1))
    return a 

def success_prob_3d(PSI):
    prob = 0
    for i in range(PSI.shape[0]): 
        if i%2 != 0 :
            prob = prob + abs(PSI[i])**2
    # prob = abs(PSI[1])**2 + abs(PSI[3])**2 
    return prob[0][0]




import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get3d_data():
    ggm_3d = np.zeros((50,50))
    for i in range(50):
        for j in range(50):
            ggm_ij = HHL_negativity.New_GGM(PSI_3d_GGM(np.sqrt(i/100.0),np.sqrt(j/100.0)))
            ggm_3d[i,j] = ggm_ij
            print(np.sqrt(i/100.0),np.sqrt(j/100.0),ggm_ij)
    return ggm_3d

def write3d_data():
    ggm_3d = np.zeros((50, 50))
    sp_3d = np.zeros((50, 50))

    with open("hhl_3d_data_sp.txt", "w") as file:
        for i in range(50):
            for j in range(50):
                ggm_ij = HHL_negativity.New_GGM(PSI_3d_GGM(np.sqrt(i/100.0), np.sqrt(j/100.0)))
                ggm_3d[i, j] = ggm_ij
                sp_ij = success_prob_3d(PSI_3_3d_GGM(np.sqrt(i/100.0), np.sqrt(j/100.0)))
                sp_3d[i,j] = sp_ij
                file.write(f"{i} {j} {ggm_ij} {sp_ij}\n")
            print(i)
        file.close()
    return ggm_3d,sp_3d

def read3d_data():
    ggm_3d = np.zeros((50, 50))
    with open("hhl_3d_data.txt", "r") as file:
        lines = file.readlines()
        for line in lines:
            values = line.strip().split()
            i = int(float(values[0]))
            j = int(float(values[1]))
            ggm_ij = float(values[2])
            ggm_3d[i, j] = ggm_ij
    return ggm_3d

def plot_3d(ggm):
    ggm_np = ggm
    plt.rcParams["text.usetex"] = True
    # Create a meshgrid for x and y coordinates
    x = np.arange(0, ggm_np.shape[1])
    y = np.arange(0, ggm_np.shape[0])
    
    # x_sq = (x/100.0)**2
    # y_sq = (y/100.0)**2

    x_sq = (x/100.0)
    y_sq = (y/100.0)


    X, Y = np.meshgrid(x_sq, y_sq)

    # Create a figure and axis object
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(X, Y, ggm_np, cmap='viridis')

    # Set labels and title
    ax.set_xlabel(r"$b_{0}^{2}$",  fontsize=16)
    ax.set_ylabel(r"$b_{1}^{2}$",  fontsize=16)

    ax.set_zlabel(r'$\mathcal{E}$', fontsize=16)
    ax.set_title(r'3D GGM')

    # Show the plot
    plt.show()
    return
