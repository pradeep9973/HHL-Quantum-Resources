# all the functions at one place

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
import hhl_3d
import hhl_coherence
import hhl_consistent
import HHL_negativity 
import hhl_4d_test as hhl4d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#success_prob in hhl_3d

# for 4d state use hhl4d.PSI_3_4d_GGM

def get_psi2_2d(b_0,eta_1 =1, eta_2 =2):
    print("*****************", b_0**2)
    return HHL_negativity.PSI_GGM(b_0,eta_1,eta_2)

def get_psi3_2d(b_0,eta_1 =1, eta_2 =2):
    return HHL_negativity.PSI_3(b_0,eta_1,eta_2)

def get_psi2_4d(b_0,b_1,b_2,eta_1 =1, eta_2 =2,eta_3 = 3,eta_4 = 4):
    return hhl4d.PSI2_4d_GGM(b_0,b_1,b_2,eta_1, eta_2 ,eta_3 ,eta_4 )

def get_psi3_4d(b_0,b_1,b_2,eta_1 =1, eta_2 =2,eta_3 = 3,eta_4 = 4):
    return hhl4d.PSI3_4d_GGM(b_0,b_1,b_2,eta_1 , eta_2 ,eta_3 ,eta_4 )

def get_psi2_3d(b_0,b_1,eta_1 =1, eta_2 =2,eta_3 = 3):
    return hhl_3d.PSI_3d_GGM(b_0,b_1,eta_1, eta_2,eta_3)

def get_psi3_3d(b_0,b_1,eta_1 =1, eta_2 =2,eta_3 = 3):
    return hhl_3d.PSI_3_3d_GGM(b_0,b_1,eta_1, eta_2,eta_3)

def get_neg(PSI, sys , sub_sys ):
    # sys should be the index which remains after partial trace
    # sub_sys is the list which indicate which should be transposed 
    # ex. [1,0,1] this means first and last index should be transposed
    return hhl_consistent.custom_negativity(PSI, sys , sub_sys )

def get_lneg(PSI, sys, sub_sys):
    return np.log2( 2*get_neg(PSI, sys,sub_sys) + 1)

def get_sp(psi_3):
    #input the PSI_3
    return hhl_3d.success_prob_nd(psi_3)

def get_coherence(psi_dm):
    # psi_dm = qutip.ket2dm(psi)
    return hhl_coherence.get_coherence(psi_dm)

def get_GGM(PSI):
    return HHL_negativity.New_GGM(PSI)

def get_GGM_AU_R(PSI):
    return


def plot_3d_sp():
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
            # psi = hhl_consistent.PSI_3d_GGM_2Din3D(b_0, b_1, 1,2,3)
            psi = get_psi2_3d(b_0,b_1)
            print(X[j, i],Y[j, i])
            # Z[j, i] = HHL_negativity.New_GGM(psi)
            Z[j, i] = get_coherence(qutip.ket2dm(psi).ptrace([2]))
            

    # Create the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')

    # Set labels and title
    ax.set_xlabel('b_0^2')
    ax.set_ylabel('b_1^2')
    ax.set_zlabel('C_r')
    ax.set_title('3D Plot of my_function')

    # Show the plot
    plt.show()

def tr_rho_sq(C,b,l1,l2):
    t1 = (1-C**2/l2**2 + b**2*(C**2/l2**2 - C**2/l1**2))**2
    # print(t1)
    t2 = (np.sqrt(1-C**2/l2**2)*C/l2 + b**2*( np.sqrt(1-C**2/l1**2)*C/l1 - np.sqrt(1-C**2/l2**2)*C/l2  )  )**2
    # print(t2)
    t3 = (C**2/l2**2 + b**2*(C**2/l1**2 - C**2/l2**2))**2
    # print(t3)
    return t1+2*t2+t3


def plot_3d_data(data_3d,title = r"$\mathcal{SP}$"):


    # ggm_np = data_3d
    plt.rcParams["text.usetex"] = True
    # Create a meshgrid for x and y coordinates
    x = np.arange(0, data_3d.shape[1])
    y = np.arange(0, data_3d.shape[0])
    
    # x_sq = (x/100.0)**2
    # y_sq = (y/100.0)**2

    x_sq = (x/100.0)
    y_sq = (y/100.0)


    X, Y = np.meshgrid(x_sq, y_sq)
    Z = data_3d
    # Create a figure and axis object
    plt.contourf(X, Y, Z,  cmap='viridis')

    # Set labels and title
    plt.xlabel(r"$b_{0}^{2}$")
    plt.ylabel(r"$b_{1}^{2}$")
    plt.title(title)

    # Add a color bar for reference
    plt.colorbar()

    # Show the plot
    plt.show()
    return

def get_3d_data():
    data_3d = np.zeros((50, 50))
    for i in range(50):
        for j in range(50):
            # ggm_ij = HHL_negativity.New_GGM(PSI_3_3d_GGM_checking(np.sqrt(i/100.0), np.sqrt(j/100.0)))
            # ggm_3d[i, j] = ggm_ij
            # sp_ij = success_prob_3d(PSI_3_3d_GGM(np.sqrt(i/100.0), np.sqrt(j/100.0)))
            sp_ij = get_coherence(qutip.ket2dm(get_psi2_3d(np.sqrt(i/100.0), np.sqrt(j/100.0))))/23.0
            # sp_ij = get_lneg(get_psi2_3d(np.sqrt(i/100.0), np.sqrt(j/100.0)),[0,1],[0,1] )
            
            data_3d[i,j] = sp_ij
            # file.write(f"{i} {j} {ggm_ij} {sp_ij}\n")
        print(i)
    return data_3d

def write_3d_data(file_name):
    data = get_3d_data()
    with open(file_name,"w") as file : 
            for i in range(50):
                for j in range(50):
                    file.write(f"{i} {j} {data[i,j]}\n")
                # print(i)
    return

def write_2d_data():
    coh_r = [get_coherence(get_psi2_2d(np.sqrt(i/100.0)).ptrace([2])) for i in range(101)]
    coh_full = [get_coherence(qutip.ket2dm(get_psi2_2d(np.sqrt(i/100.0))))/15.0 for i in range(101)]
    ggm2d = [get_GGM(get_psi2_2d(np.sqrt(i/100.0))) for i in range(101) ]
    lneg2d = [get_lneg(get_psi2_2d(np.sqrt(i/100.0)),[0,1],[0,1]) for i in range(101)]
    b0_sq = [i/100.0 for i in range(101)]

    with open("Data_files/2D/ggm2d.txt","w") as file:
        for i in ggm2d:
            file.write(str(i) + "\n")

    with open("Data_files/2D/lneg2d.txt","w") as file:
        for i in lneg2d:
            file.write(str(i) + "\n")    

    with open("Data_files/2D/C_r_2d.txt","w") as file:
        for i in coh_r:
            file.write(str(i) + "\n")       

    with open("Data_files/2D/C_full_2d.txt","w") as file:
        for i in coh_full:
            file.write(str(i) + "\n")

    with open("Data_files/2D/b0_sq.txt","w") as file:
        for i in b0_sq:
            file.write(str(i) + "\n")