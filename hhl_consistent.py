# checking consistency


import numpy as np
from numpy.linalg import eig

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
import hhl_3d

# Define a 6x6 matrix
# A = np.array([[14/6.0, -4/6.0, -4/6.0, 0, 0, 0],
#               [-4/6.0, 11/6.0, -1/6.0, 0, 0, 0],
#               [-4/6.0, -1/6.0, 11/6.0, 0, 0, 0],
#               [ 0, 0, 0, 14/6.0, -4/6.0, -4/6.0],
#               [ 0, 0, 0,-4/6.0, 11/6.0, -1/6.0],
#               [ 0, 0, 0,-4/6.0, -1/6.0, 11/6.0]])





# # Find the eigenvalues and eigenvectors of A
# eigvals, eigvecs = np.linalg.eigh(A)




def PSI_3_3d_GGM(b_0=0.5,b_1=0.4,eta_1=1,eta_2=2,eta_3=3):
    # here the eigenvalue is represented by a single vector

    
    b_2 = np.sqrt(1 - b_0**2 - b_1**2)
    print(b_2)



    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    beta_1 = (b_0 - b_1 )/np.sqrt(2)
    beta_2 = (b_0 + b_1)/np.sqrt(2)
    beta_3 = (b_2)/np.sqrt(1)

    

    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736

    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(3,0) - basis(3,1) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(3,0) + basis(3,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (basis(3,2)).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    PSI = beta_1*tensor(psi_c_1,psi_d_1) + beta_2*tensor(psi_c_2,psi_d_2) + beta_3*tensor(psi_c_3,psi_d_3)
    # print(PSI)
    return PSI


def success_prob_3d(PSI):
    prob = 0
    for i in range(PSI.shape[0]): 
        if i%2 != 0 :
            prob = prob + abs(PSI[i])**2
    # prob = abs(PSI[1])**2 + abs(PSI[3])**2 
    return prob[0][0]

def PSI_3d_GGM_2Din3D(b_0=0.2,b_1=0.4,eta_1=1,eta_2=2,eta_3=3):
    # here the eigenvalue is represented by a single vector

    
    b_2 = np.sqrt(1 - b_0**2 - b_1**2)
    # print(b_2)


    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    beta_1 = (b_0 + b_1)/np.sqrt(2)
    beta_2 = (b_0 - b_1)/np.sqrt(2)
    beta_3 = (b_2)/np.sqrt(1)

    

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

    psi_c_1 = (basis(3,0) + basis(3,1)).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(3,0) - basis(3,1)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (basis(3,2)).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    PSI = beta_1*tensor(psi_eta_1,psi_c_1,psi_d_1) + beta_2*tensor(psi_eta_2,psi_c_2,psi_d_2) + beta_3*tensor(psi_eta_3,psi_c_3,psi_d_3)
    # print(PSI)
    return PSI






import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_3d():
    # Generate data points
    x = np.linspace(0, 0.49, 25)
    y = np.linspace(0, 0.49, 25)
    X, Y = np.meshgrid(x, y)

    # Create an empty Z array
    Z = np.zeros_like(X)

    # Populate Z array using a loop
    for i in range(len(x)):
        for j in range(len(y)):
            b_0 = np.sqrt(X[j, i])
            b_1 = np.sqrt(Y[j, i])
            b_2 = np.sqrt(1-b_0**2-b_1**2)
            psi = PSI_3d_GGM_2Din3D(b_0, b_1, 1,2,3)
            # psi = hhl_3d.PSI_3d_GGM(b_0,b_1)
            print(X[j, i],Y[j, i])
            Z[j, i] = HHL_negativity.New_GGM(psi)

    # Create the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')

    # Set labels and title
    ax.set_xlabel('b_0^2')
    ax.set_ylabel('b_1^2')
    ax.set_zlabel('sp')
    ax.set_title('3D Plot of my_function')

    # Show the plot
    plt.show()

def custom_negativity(PSI, sys , sub_sys ):
    if len(sys) == 3:
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


def plot_3d_negativity(sub_sys,sys_ptrace):
    # Generate data points
    x = np.linspace(0, 0.49, 25)
    y = np.linspace(0, 0.49, 25)
    X, Y = np.meshgrid(x, y)

    # Create an empty Z array
    Z = np.zeros_like(X)

    # Populate Z array using a loop
    for i in range(len(x)):
        for j in range(len(y)):
            b_0 = np.sqrt(X[j, i])
            b_1 = np.sqrt(Y[j, i])
            b_2 = np.sqrt(1-b_0**2-b_1**2)
            # psi = PSI_3d_GGM_2Din3D(b_0, b_1, 1,2,1 )
            psi = hhl_3d.PSI_3d_GGM(b_0,b_1)
            print(X[j, i],Y[j, i])
            Z[j, i] = custom_negativity(psi,sub_sys,sys_ptrace)

    # Create the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')

    # Set labels and title
    ax.set_xlabel('b_0^2')
    ax.set_ylabel('b_1^2')
    ax.set_zlabel('Log_N')
    ax.set_title('3D Plot of my_function')

    # Show the plot
    plt.show()
