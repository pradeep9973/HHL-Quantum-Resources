#author: pradeep kumar
# module : GGM

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


#schmidt decompostion

def get_binary(n,pad=0):
    b = bin(n)[2:]
    if pad > len(b):
        b = "0"*(pad - len(b)) + b
    return b
        


# def getting_binary_list(length):
#     digit_length = len(bin(length)) - 2
#     key = [ "0"*(digit_length - int())  bin(i)[2:] for i in range(length)]

def single_states_matrix(PSI, k):
    # n_qubit_state to schmidt matrix of kth qubit
    PSI = np.array(PSI)
    key = [get_binary(i,int(np.log2(len(PSI)))) for i in range(len(PSI))]
    # print(key)
    schmidt_mat = np.zeros([2,int(len(PSI)/2)],dtype=complex)
    # print(schmidt_mat)
    i_0 = 0
    i_1 = 0 
    for i in range(len(key)):
        # print(key[i],PSI[i])
        # schmidt_mat[0][i_0] = PSI[i]
        if int(key[i][k]) == 0 :
            schmidt_mat[0][i_0] = PSI[i]
            i_0 = i_0 + 1
        else :
            schmidt_mat[1][i_1] = PSI[i]
            i_1 = i_1 + 1
    return schmidt_mat

# def two_state_matrix()

def all_states_matrix(PSI, index_list = [0]):
    # this generates the schimdt matrix of the state, whose eigenvalues gives the lambda's of the schimdt decompositon
    PSI = np.array(PSI)
    key = [get_binary(i,int(np.log2(len(PSI)))) for i in range(len(PSI))]

    schmidt_mat = np.zeros([2**len(index_list),int(len(PSI)/(2**len(index_list)))],dtype=complex)
    # print("ERROR")
    # print(schmidt_mat)
    i_0 = 0
    i_1 = 0 
    matrix_rows = [get_binary(i,len(index_list)) for i in range(2**len(index_list))]
    for i in range(len(matrix_rows)) :
        dummy = 0
     
        for k in range(len(key)):
            sub_key = [key[k][i] for i in index_list]
            if matrix_rows[i] == "".join(sub_key):
                schmidt_mat[i][dummy] = PSI[k]
                dummy = dummy + 1

    return schmidt_mat



