#author: pradeep kumar
# plotting two curves with different set of eigenvalues


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
import HHL_negativity
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

def GGM_plots(eigs_1 = [1,2],eigs_2 = [3,5]):
    lambda_1s = HHL_negativity.get_GGM(eigs_1[0],eigs_1[1])
    lambda_2s = HHL_negativity.get_GGM(eigs_2[0],eigs_2[1])
    
    plt.rcParams["text.usetex"] = True

    fig, ax = plt.subplots()

    # ax.plot(lambda_1_2[0], lambda_1_2[1], label=r'$\mathrm{\frac{\lambda_{1}}{\lambda_{2}} = \frac{1}{2}}$')
    ax.plot(lambda_1s[0], lambda_1s[1], label=r'$\mathrm{\frac{\lambda_{1}}{\lambda_{2}} = \frac{' + str(eigs_1[0]) +'}{'+ str(eigs_1[1])+ '}}$')

    # ax.plot(lambda_2_5[0], lambda_2_5[1], label=r'$\mathrm{\frac{\lambda_{1}}{\lambda_{2}} = \frac{2}{5}}$')
    ax.plot(lambda_2s[0], lambda_2s[1], label=r'$\mathrm{\frac{\lambda_{1}}{\lambda_{2}} = \frac{' + str(eigs_2[0]) +'}{'+ str(eigs_2[1])+ '}}$')

    
   
    ax.grid(True)

    ax.legend()
    ax.set_xlabel(r'$\mathrm{b}_{\mathrm{0}}^{\mathrm{2}}$')

    ax.set_ylabel(r'GGM')
    plt.title(r"GGM vs $\mathrm{b}_{\mathrm{0}}^{\mathrm{2}}$")
    # plt.title(r"GGM vs $\mathrm{\lambda_{1}$")
    # plt.title(r'GGM vs $\mathrm{\lambda_{1}}$')

    plt.savefig("figure.pgf")
    plt.show()


def plot_negativity():
    eta_1 = 7
    eta_2 = 8
    lambda_1_2 = HHL_negativity.get_negativity(eta_1,eta_2)
    lambda_2_5 = HHL_negativity.get_negativity(2,3)
    lambda_4_5 = HHL_negativity.get_negativity(4,5)
    lambda_2_4 = HHL_negativity.get_negativity(2,4)

    # lambda_2_5 = HHL_negativity.get_GGM(2,5)

    plt.rcParams["text.usetex"] = True

    fig, ax = plt.subplots()


    ax.plot(lambda_1_2[0], lambda_1_2[1], label=r'$\mathrm{\frac{\lambda_{2}}{\lambda_{1}} = \frac{2}{1}}$')
    ax.plot(lambda_2_5[0], lambda_2_5[1], label=r'$\mathrm{\frac{\lambda_{2}}{\lambda_{1}} = \frac{3}{2}}$')
    ax.plot(lambda_4_5[0], lambda_4_5[1], label=r'$\mathrm{\frac{\lambda_{2}}{\lambda_{1}} = \frac{5}{4}}$')
    ax.plot(lambda_2_4[0], lambda_2_4[1], label=r'$\mathrm{\frac{\lambda_{2}}{\lambda_{1}} = \frac{4}{2}}$')

    # ax.plot(lambda_1_2[0], lambda_1_2[1], label=r'$\mathrm{\lambda_{1} = 1, \lambda_{2} = 2}$')
    # ax.plot(lambda_2_5[0], lambda_2_5[1], label=r'$\mathrm{\lambda_{1} = 2, \lambda_{2} = 5}$')



    ax.grid(True)

    ax.legend()
    ax.set_xlabel(r'$\mathrm{b}_{\mathrm{0}}^{\mathrm{2}}$')

    ax.set_ylabel(r'Negativity')
    plt.title(r"$\mathrm{Negativity([\lambda][U]|AUX)}$ vs $\mathrm{b}_{\mathrm{0}}^{\mathrm{2}}$")
    # plt.title(r"GGM vs $\mathrm{\lambda_{1}$")
    # plt.title(r'GGM vs $\mathrm{\lambda_{1}}$')

    plt.savefig("figure_neg.pgf")
    plt.show()


# def plot_3d(b_0_sq, negativity, eta_1_values):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_xlabel('b_0^2')
#     ax.set_ylabel('Negativity')
#     ax.set_zlabel('Ratio of eta_1 to eta_2')
#     for i in range(len(eta_1_values)):
#         eta_1 = eta_1_values[i]
#         eta_2 = eta_1 + 1
#         x, y = HHL_negativity.get_negativity(eta_1, eta_2)
#         print("++++++++++done++++++++++",i)
#         ax.scatter(x, y, eta_1/eta_2, c='red', marker='o')
#     plt.show()

# eta_1_values = [i for i in range(1, 11)]
# plot_3d(*HHL_negativity.get_negativity(), eta_1_values)