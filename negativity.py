# author : pradeep jangra
# negativity of entanglement



import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

def partial_transpose(mat,a):
    #ensure mat is well defined
    mat_copy = np.copy(mat)
    s = int(np.sqrt(mat.shape[0]))
    for w in range(0,s):
        for v in range(0,s):
            for j in range(0,s):
                for i in range(0,s):
                    if (a == 0):
                        mat_copy[v*s + i,w*s+j] = mat[w*s + i,v*s+j]
                    elif(a==1) :
                        mat_copy[v*s + i,w*s+j] = mat[v*s + j,w*s+i]
                    # print((v*2 + i,w*2+j),mat_copy[v*2 + i,w*2+j])
    
    return mat_copy

def partial_transpose_1(mat,a,d_1,d_2):
    #d_1 dimension of first index
    #d_2 dimension of second index
    #ensure mat is well defined
    mat_copy = np.copy(mat)
    # s = int(np.sqrt(mat.shape[0]))
    for w in range(0,d_1):
        for v in range(0,d_2):
            for j in range(0,d_1):
                for i in range(0,d_2):
                    if (a == 0):
                        mat_copy[v*d_1 + i,w*d_2+j] = mat[w*d_2 + i,v*d_1+j]
                    elif(a==1) :
                        mat_copy[v*d_1 + i,w*d_2+j] = mat[v*d_1 + j,w*d_2+i]
                    # print((v*2 + i,w*2+j),mat_copy[v*2 + i,w*2+j])
    
    return mat_copy

def negativity(rho,index_no = 0):
    rho_p_t = partial_transpose(rho,index_no)
    return (LA.norm(rho_p_t, "nuc") - 1)/2.0

def negativity_direct(rho):

    return (LA.norm(rho, "nuc") - 1)/2.0



test_rho = np.loadtxt("rho_ab.txt",dtype=float)
# print(negativity(rho_ab,0))
# print(negativity(rho_ab,1))

# calulating N_ab
b_0 = 0.2
b_1 = np.sqrt(1-b_0**2)

z = np.array([[1],[0]])
i = np.array([[0],[1]])


eta_1 = 1
eta_2 = 2
C = 0.736

b_m = (b_0-b_1)**2
b_p = (b_0+b_1)**2

# rho_ab = p_1Xaux_1

p_1 = z - i
p_2 = z + i
aux_1 = np.sqrt(1-C**2/eta_1**2)*z + (C/eta_1)*i
aux_2 = np.sqrt(1-C**2/eta_2**2)*z + (C/eta_2)*i

# print(aux_1)

p_aux_1 = np.kron(p_1, aux_1)
p_aux_2 = np.kron(p_2,aux_2)


rho_ab = (b_m/4.0)* (np.kron(np.kron(p_1,aux_1),np.kron(p_1,aux_1).T)) + b_p/4.0*(np.kron(np.kron(p_2,aux_2),np.kron(p_2,aux_2).T))

q_1 = np.array([[0],[1],[0],[0]])
q_2 = np.array([[0],[0],[1],[0]])


rho_ac = (b_m/2.0 )* (np.kron(np.kron(q_1,aux_1),np.kron(q_1,aux_1).T)) + (b_p/2.0)* (np.kron(np.kron(q_2,aux_2),np.kron(q_2,aux_2).T))

print(rho_ab)
pt_ab = partial_transpose(rho_ab,0)
print(pt_ab)
# print(rho_ac)

print("Negativity is : ",negativity(rho_ab))
print(LA.eigh(pt_ab))
print(negativity_direct(rho_ac))

# print((LA.eigh(test_rho)[0]))
# print(negativity(test_rho))
# print("$$")
# print(LA.norm(test_rho, "nuc"))
# print(np.sum(np.linalg.svd(test_rho)[1]))

print("===========================")
# calculation N_a_bc

# three_010 = np.array([[0],[0],[1],[0],[0],[0],[0],[0]])
# three_011 = np.array([[0],[0],[0],[1],[0],[0],[0],[0]])
# three_100 = np.array([[0],[0],[0],[0],[1],[0],[0],[0]])
# three_101 = np.array([[0],[0],[0],[0],[0],[1],[0],[0]])


# psi_abc = ((b_0-b_1)/2.0) * np.kron((three_010-three_011),aux_1) + ((b_0+b_1)/2.0) * np.kron((three_100 - three_101),aux_2)
# rho_a_bc = np.kron(psi_abc,psi_abc.T)
# print(negativity_direct(rho_a_bc))