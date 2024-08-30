# calculate and plot coherece for n-D states

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
import hhl_consistent
import hhl_4d_test





def get_coherence(psi_dm):

    # Convert the pure state to a density matrix
    # density_matrix = qutip.ket2dm(psi)
    density_matrix = psi_dm
    # return density_matrix
    shape = density_matrix.shape
    rows = shape[0]
    columns = shape[1]

    # Calculate the coherence measure (L1 norm)
    l1_norm = 0
    for i in range(rows):
        for j in range(columns):
            if i != j :
                l1_norm += np.abs(density_matrix[i,j])
        
    coherence_measure = l1_norm

    # Print the coherence measure
    # print("L1 norm of coherence:", coherence_measure)
    return coherence_measure

def coherence_2d(b_0):
    psi = HHL_negativity.PSI_GGM(b_0)
    coherence = get_coherence(qutip.ket2dm(psi).ptrace([2]))
    return coherence

def coherence_2d_local(b_0, sub_system):
    psi = HHL_negativity.PSI_GGM(b_0)
    psi_dm = qutip.ket2dm(psi)
    coherence = get_coherence(psi_dm.ptrace([sub_system]))
    return coherence


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_3d_coherence(sub_sys):

    # Generate data points
    x = np.linspace(0, 0.49, 50)
    y = np.linspace(0, 0.49, 50)
    X, Y = np.meshgrid(x, y)

    # Create an empty Z array
    Z = np.zeros_like(X)

    # Populate Z array using a loop
    for i in range(len(x)):
        for j in range(len(y)):
            b_0 = np.sqrt(X[j, i])
            b_1 = np.sqrt(Y[j, i])
            # b_2 = np.sqrt(1-b_0**2-b_1**2-0.1)
            # b_2 = 0
            # psi = hhl_3d.PSI_3_3d_GGM_checking(b_0, b_1)
            psi = hhl_consistent.PSI_3d_GGM_2Din3D(b_0, b_1,1,2,3)

            psi_dm = qutip.ket2dm(psi)
            Z[j, i] = get_coherence(psi_dm.ptrace(sub_sys))

    # Create the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')

    # Set labels and title
    ax.set_xlabel('b_0^2')
    ax.set_ylabel('b_1^2')
    ax.set_zlabel('coherence')
    ax.set_title('Coherence vs b_0^2')

    # Show the plot
    plt.show()


def plot_2d_contour_coherence(sub_sys):

    # Generate data points
    x = np.linspace(0, 0.49, 50)
    y = np.linspace(0, 0.49, 50)
    X, Y = np.meshgrid(x, y)

    # Create an empty Z array
    Z = np.zeros_like(X)

    # Populate Z array using a loop
    for i in range(len(x)):
        for j in range(len(y)):
            b_0 = np.sqrt(X[j, i])
            b_1 = np.sqrt(Y[j, i])
            # b_2 = np.sqrt(1-b_0**2-b_1**2-0.1)
            # b_2 = 0
            psi = hhl_3d.PSI_3d_GGM(b_0, b_1)
            # psi = hhl_consistent.PSI_3d_GGM_2Din3D(b_0, b_1,1,2,3)

            psi_dm = qutip.ket2dm(psi)
            Z[j, i] = get_coherence(psi_dm.ptrace(sub_sys))

    # Create the 3D plot
    # ll = 4
    # ul = 11
    # contour_line_values = [ (ll + (ul - ll)*i/100.0) for i in range(100) ]

    # plt.contourf(X, Y, Z, contour_line_values, cmap='viridis')
    plt.contourf(X, Y, Z, cmap='viridis')
    

    # Set labels and title
    plt.xlabel('b_0^2')
    plt.ylabel('b_1^2')
    plt.title('Coherence vs b_0^2')

    # Add a color bar for reference
    plt.colorbar()

    # Show the plot
    plt.show()


