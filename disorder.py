# author : pradeep
# disorder in HHL

import random
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
import HHL_negativity

# Generate a random number from a Gaussian distribution with mean 10 and standard deviation 2
random_number = random.gauss(0, 0.02)

# print(random_number)



def PSI_3(b_0,eta_1=1,eta_2=2):
    print("b_0 : ", b_0)
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

def PSI_3_disorder(b_0,eta_1=1,eta_2=2,w=0.1,u=0.2):
    # print(b_0)


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

def x_disorder (b0,l1,l2,w):
    b1 = np.sqrt(1-b0**2)
    xa = 0.5*( ((b0 - b1)/l1)*np.exp(complex(0,w)) + (b0 + b1)/l2 )
    xb = 0.5*( (-1*(b0 - b1)/l1)*np.exp(complex(0,w)) + (b0 + b1)/l2 )
    arr = np.array([xa,xb])
    # print(arr)
    return arr

def x_original(b0, l1, l2):
    return x_disorder(b0,l1, l2, 0)

def succ_distance(psi3_og, psi3_d):

    # return abs(np.dot(a1.conj(),a2)/(abs(a1)*abs(a2)))
    # return abs(np.dot(a1.conj(),a2) - np.dot(a2.conj(),a2))
    a1 = np.array([psi3_og[1][0][0],psi3_og[3][0][0]])
    a2 = np.array([psi3_d[1][0][0],psi3_d[3][0][0]])

    # print(np.linalg.norm(a1),np.linalg.norm(a2))
    return np.linalg.norm(a1-a2)/np.linalg.norm(a1)

def state_fidelity(psi3_og, psi3_d):
    # return abs(np.dot(a1.conj(),a2)/(abs(a1)*abs(a2)))
    # return abs(np.dot(a1.conj(),a2) - np.dot(a2.conj(),a2))
    a1 = np.array([psi3_og[1][0][0],psi3_og[3][0][0]])
    a2 = np.array([psi3_d[1][0][0],psi3_d[3][0][0]])

    # print(np.linalg.norm(a1),np.linalg.norm(a2))
    return np.linalg.norm(np.dot(np.conj(a1), a2))


def success_prob(PSI):
    prob = abs(PSI[1])**2 + abs(PSI[3])**2 
    return prob[0][0]

# def plot_efficiency_vs_b_0_sq():
#     m = []
#     B  = [i/100.0 for i in range(0,101)]
#     B_sq = [b**2 for b in B]
#     actual_prob = []
#     for b in B :
#         d = []
#         for i in range(1): 
#             w = random.gauss(0, 3)
#             # dis = distance(x_disorder(b,1,2,w), x_original(b,1,2))
#             dis = success_prob(PSI_3_disorder(b,1,2,w))
#             # print(i, " ====== " , dis)
#             d.append(dis)
#         avg = np.mean(d)
#         m.append(avg)
#         actual = success_prob(PSI_3(b,1,2))
#         actual_prob.append(actual)
#         print(b,"-->", actual,avg)

#     # for b in B :
#     #     d = []
#     #     for i in range(100): 
#     #         w = random.gauss(0, 0.1)
#     #         dis = distance(x_disorder(b,1,2,w), x_original(b,1,2))
#     #         # print(i, " ====== " , dis)
#     #         d.append(dis)
#     #     avg = np.mean(d)
#     #     m.append(avg)
#     #     print(b,"-->", avg)

#     # plt.plot(B_sq,m,label="disorder")

#     plt.rcParams["text.usetex"] = True
#     fig, ax = plt.subplots()
#     x = B_sq
#     y1 = actual_prob
#     y2 = m
#     ax.plot(x, y1, label=r'Efficiency',color='#ff6700')
#     ax.plot(x, y2, label=r'disorder')

#     # Add grid and legend
#     # ax.grid(True)
#     ax.legend()

#     # Set title and axes labels in LaTeX
#     # ax.set_title(r'\textbf{Entanglement vs ratio of eigenvalues}')
#     # ax.set_xlabel(r"$\mathrm{\frac{\lambda_{2}}{\lambda_{1}}}$")
#     ax.set_xlabel(r"$b_{0}^{2}$")

#     ax.set_ylabel(r'Efficiency')
#     # fig.savefig("fig_3.eps", format="eps")
#     plt.show()
#     return

# def efficiency_vs_ratio():
#     # m = []
#     b = 0.1
#     ratio = [i/100.0 for i in range(1,100)]
#     # B_sq = [b**2 for b in B]
#     actual_prob = []
#     for r in ratio :
#         d = []
#         # for i in range(100): 
#         #     w = random.gauss(0, 0.3)
#         #     # dis = distance(x_disorder(b,1,2,w), x_original(b,1,2))
#         #     dis = success_prob(PSI_3_disorder(b,1,2,w))
#         #     # print(i, " ====== " , dis)
#         #     d.append(dis)
#         # avg = np.mean(d)
#         # m.append(avg)
#         actual = success_prob(PSI_3(b,int(r*100),100))
#         actual_prob.append(actual)
#         print(b,"-->", actual)

#     # for b in B :
#     #     d = []
#     #     for i in range(100): 
#     #         w = random.gauss(0, 0.1)
#     #         dis = distance(x_disorder(b,1,2,w), x_original(b,1,2))
#     #         # print(i, " ====== " , dis)
#     #         d.append(dis)
#     #     avg = np.mean(d)
#     #     m.append(avg)
#     #     print(b,"-->", avg)

#     # plt.plot(B_sq,m,label="disorder")


#     psi_list = [HHL_negativity.PSI_GGM(0.2,int(i*100),100) for i in ratio]
    
#     ggm =[]
#     counter = 0
#     for i in psi_list:
#         g = HHL_negativity.New_GGM(i)
#         ggm.append(g)
#         print(counter,g)
#         counter+=1
#     plt.rcParams["text.usetex"] = True
#     fig, ax = plt.subplots()
#     x = ratio
#     y1 = actual_prob
#     y2 = ggm
#     ax.plot(x, y1, label=r'Efficiency',color='black')
#     ax2 = ax.twinx()
#     ax2.plot(x, y2, label=r'GGM', color="black", linestyle='dashed')

#     # Add grid and legend
#     # ax.grid(True)
#     # ax.legend()

#     # Set title and axes labels in LaTeX
#     # ax.set_title(r'\textbf{Entanglement vs ratio of eigenvalues}')
#     ax.set_xlabel(r"$\mathrm{\frac{\lambda_{2}}{\lambda_{1}}}$")
#     # ax.set_xlabel(r"$b_{0}^{2}$")

#     ax.set_ylabel(r'Efficiency')
#     ax2.set_ylabel(r'GGM')

#     # fig.savefig("fig_3.eps", format="eps")
#     plt.show()
