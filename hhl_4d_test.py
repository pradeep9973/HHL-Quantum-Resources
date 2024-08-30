# testing 4 dimensions

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


def PSI2_4d_GGM(b_0=0.5,b_1=0.2,b_2=0.1,eta_1=1,eta_2=2,eta_3=3, eta_4=4):
    # here the eigenvalue is represented by a single vector

    
    b_3 = np.sqrt(1 - b_0**2 - b_1**2 - b_2**2)
    # print(b_2)


    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    beta_1 = (-b_0 + b_1 + b_2 + b_3)/np.sqrt(4)
    beta_2 = (b_0 - b_1 + b_2 + b_3)/np.sqrt(4)
    beta_3 = (b_0 + b_1 - b_2 + b_3)/np.sqrt(4)
    beta_4 = (b_0 + b_1 + b_2 - b_3)/np.sqrt(4)

    

    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736
    eta1_bin = bin(eta_1)[2:]
    eta2_bin = bin(eta_2)[2:]
    eta3_bin = bin(eta_3)[2:]
    eta4_bin = bin(eta_4)[2:]


    psi_eta_1 = tensor(basis(1))
    psi_eta_2 = tensor(basis(1))
    psi_eta_3 = tensor(basis(1))
    psi_eta_4 = tensor(basis(1))

    len_lambda = max(len(eta1_bin),len(eta2_bin),len(eta3_bin),len(eta4_bin))

 
    if len(eta1_bin) < len_lambda:
        eta1_bin = "0"*(len_lambda - len(eta1_bin)) + eta1_bin
    if len(eta2_bin) < len_lambda:
        eta2_bin = "0"*(len_lambda - len(eta2_bin)) + eta2_bin
    if len(eta3_bin) < len_lambda:
        eta3_bin = "0"*(len_lambda - len(eta3_bin)) + eta3_bin
    if len(eta4_bin) < len_lambda:
        eta4_bin = "0"*(len_lambda - len(eta4_bin)) + eta4_bin
    # print(eta1_bin,eta2_bin)  
    # print(eta1_bin,eta2_bin,eta3_bin)  
 
    # len_lamda = len(eta2_bin)

    for i in range(len_lambda):

        psi_eta_1 = tensor(basis(2,int(eta1_bin[len_lambda-i-1])),psi_eta_1)
        psi_eta_2 = tensor(basis(2,int(eta2_bin[len_lambda-i-1])),psi_eta_2)
        psi_eta_3 = tensor(basis(2,int(eta3_bin[len_lambda-i-1])),psi_eta_3)
        psi_eta_4 = tensor(basis(2,int(eta4_bin[len_lambda-i-1])),psi_eta_4)


    # print(psi_eta_1,psi_eta_2)
    psi_eta_1 = Qobj(psi_eta_1.full().flatten().reshape(psi_eta_1.shape))
    psi_eta_2 = Qobj(psi_eta_2.full().flatten().reshape(psi_eta_2.shape))
    psi_eta_3 = Qobj(psi_eta_3.full().flatten().reshape(psi_eta_3.shape))
    psi_eta_4 = Qobj(psi_eta_4.full().flatten().reshape(psi_eta_4.shape))


    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (-basis(4,0) + basis(4,1) + basis(4,2) + basis(4,3) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(4,0) - basis(4,1) + basis(4,2) + basis(4,3) ).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (basis(4,0) + basis(4,1) - basis(4,2) + basis(4,3) ).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    psi_c_4 = (basis(4,0) + basis(4,1) + basis(4,2) - basis(4,3) ).unit()
    psi_d_4 = np.sqrt(1-C**2/eta_4**2)*basis(2,0) + C/eta_4*basis(2,1)

    PSI = beta_1*tensor(psi_eta_1,psi_c_1,psi_d_1) + beta_2*tensor(psi_eta_2,psi_c_2,psi_d_2) + beta_3*tensor(psi_eta_3,psi_c_3,psi_d_3) + beta_4*tensor(psi_eta_4,psi_c_4,psi_d_4)
    # print(PSI)
    return PSI




def PSI3_4d_GGM(b_0=0.5,b_1=0.4,b_2 =0.2, eta_1=1,eta_2=2,eta_3=3, eta_4 = 4):
    # here the eigenvalue is represented by a single vector

    
    b_3 = np.sqrt(1 - b_0**2 - b_1**2 - b_2**2 )
    # b_3 = np.sqrt(0.1)
    # print(b_2) "{:.3f}".format(number)
    print("{:.3f}".format(b_0**2))



    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    # beta_1 = (b_0 + b_1 + b_2)/np.sqrt(3)
    # beta_2 = (b_1 - b_2)/np.sqrt(2)
    # beta_3 = (-2*b_0 + b_1 + b_2)/np.sqrt(6)

    beta_1 = (-b_0 + b_1 + b_2 + b_3)/np.sqrt(4)
    beta_2 = (b_0 - b_1 + b_2 + b_3)/np.sqrt(4)
    beta_3 = (b_0 + b_1 - b_2 + b_3)/np.sqrt(4)
    beta_4 = (b_0 + b_1 + b_2 - b_3)/np.sqrt(4)


    C = 0.736


    psi_c_1 = (-basis(4,0) + basis(4,1) + basis(4,2) + basis(4,3) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(4,0) - basis(4,1) + basis(4,2) + basis(4,3) ).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (basis(4,0) + basis(4,1) - basis(4,2) + basis(4,3) ).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    psi_c_4 = (basis(4,0) + basis(4,1) + basis(4,2) - basis(4,3) ).unit()
    psi_d_4 = np.sqrt(1-C**2/eta_4**2)*basis(2,0) + C/eta_4*basis(2,1)
    

    PSI = beta_1*tensor(psi_c_1,psi_d_1) + beta_2*tensor(psi_c_2,psi_d_2) + beta_3*tensor(psi_c_3,psi_d_3) + beta_4*tensor(psi_c_4,psi_d_4)
    # print(PSI)
    return PSI


def psi_3_consistent(b_0=0.5,b_1=0.4,b_2 =0.2, eta_1=3,eta_2=2,eta_3=1, eta_4 = 4):
    # here the eigenvalue is represented by a single vector

    
    # b_3 = np.sqrt(1 - b_0**2 - b_1**2 - b_2**2 )
    b_3 = 0
    # print(b_2)
    print(b_0**2,b_1**2)


    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    # beta_1 = (b_0 + b_1 + b_2)/np.sqrt(3)
    # beta_2 = (b_1 - b_2)/np.sqrt(2)
    # beta_3 = (-2*b_0 + b_1 + b_2)/np.sqrt(6)

    beta_1 = (b_0 + b_1 + b_2)/np.sqrt(3)
    beta_2 = (b_1 - b_2)/np.sqrt(2)
    beta_3 = (-2*b_0 + b_1 + b_2)/np.sqrt(6)
    beta_4 = (b_3)

    

    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736

    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(4,0) + basis(4,1) + basis(4,2) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(4,1) - basis(4,2)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (-2*basis(4,0) + basis(4,1) + basis(4,2)).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    psi_c_4 = (basis(4,3) ).unit()
    psi_d_4 = np.sqrt(1-C**2/eta_4**2)*basis(2,0) + C/eta_4*basis(2,1)
    

    PSI = beta_1*tensor(psi_c_1,psi_d_1) + beta_2*tensor(psi_c_2,psi_d_2) + beta_3*tensor(psi_c_3,psi_d_3) + beta_4*tensor(psi_c_4,psi_d_4)
    # print(PSI)
    return PSI


def psi_2D_4D(b_0=0.5,b_1=0.4,b_2 =0.2, eta_1=1,eta_2=1,eta_3=5, eta_4 = 5):
    # here the eigenvalue is represented by a single vector
    #[2,1]
    #[3,4]
    
    # b_3 = np.sqrt(1 - b_0**2 - b_1**2 - b_2**2 )

    b_3 = 0
    # print(b_2)
    # print(b_0**2,b_1**2)


    # b_m = (b_0-b_1)/np.sqrt(2)
    # b_p = (b_0+b_1)/np.sqrt(2)
    # beta_1 = (b_0 + b_1 + b_2)/np.sqrt(3)
    # beta_2 = (b_1 - b_2)/np.sqrt(2)
    # beta_3 = (-2*b_0 + b_1 + b_2)/np.sqrt(6)

    beta_1 = (-3*b_0 + b_1)/np.sqrt(10)
    beta_2 = (- b_2 - b_3)/np.sqrt(2)
    beta_3 = (b_0 + b_1)/np.sqrt(2)
    beta_4 = (b_2 + 3*b_3)/np.sqrt(10)

    

    # eta_1 = 1.0
    # eta_2 = 2.0
    C = 0.736

    # print(psi_eta_1,psi_eta_2)

    psi_c_1 = (basis(4,0) + basis(4,1) + basis(4,2) ).unit()
    psi_d_1 = np.sqrt(1-C**2/eta_1**2)*basis(2,0) + C/eta_1*basis(2,1)

    psi_c_2 = (basis(4,1) - basis(4,2)).unit()
    psi_d_2 = np.sqrt(1-C**2/eta_2**2)*basis(2,0) + C/eta_2*basis(2,1)

    psi_c_3 = (-2*basis(4,0) + basis(4,1) + basis(4,2)).unit()
    psi_d_3 = np.sqrt(1-C**2/eta_3**2)*basis(2,0) + C/eta_3*basis(2,1)

    psi_c_4 = (basis(4,3) ).unit()
    psi_d_4 = np.sqrt(1-C**2/eta_4**2)*basis(2,0) + C/eta_4*basis(2,1)
    

    PSI = beta_1*tensor(psi_c_1,psi_d_1) + beta_2*tensor(psi_c_2,psi_d_2) + beta_3*tensor(psi_c_3,psi_d_3) + beta_4*tensor(psi_c_4,psi_d_4)
    # print(PSI)
    return PSI


def success_prob_4d(PSI):
    prob = 0
    for i in range(PSI.shape[0]): 
        if i%2 != 0 :
            prob = prob + abs(PSI[i])**2
    # prob = abs(PSI[1])**2 + abs(PSI[3])**2 
    return prob[0][0]


def plot_4in3D():
        # import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D



    # Generate data points
    x = np.linspace(0, 0.4, 25)
    y = np.linspace(0, 0.4, 25)
    X, Y = np.meshgrid(x, y)

    # Create an empty Z array
    Z = np.zeros_like(X)

    # Populate Z array using a loop
    for i in range(len(x)):
        for j in range(len(y)):
            b_0 = np.sqrt(X[j, i])
            b_1 = np.sqrt(Y[j, i])
            b_2 = np.sqrt(1-b_0**2-b_1**2-0.1)
            # b_2 = 0
            psi = PSI_3_4d_GGM(b_0, b_1, b_2 )
            Z[j, i] = success_prob_4d(psi)

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
