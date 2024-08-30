#author : Pradeep Kumar
#summerProject


import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d











# rho_ab = np.loadtxt("rho_ab.txt",dtype=float)
# rho_aTb = np.kron(rho_a.T,rho_b)
# rho_abT = np.kron(rho_a,rho_b.T)


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

# print(rho_ab)
# rho_aTb = partial_transpose(rho_ab,0)
# print(rho_aTb)
#0 for aT
#1 for bT

# trace norm = LA.norm(rho_aTb, "nuc")

def negativity(aTb):
    return (LA.norm(aTb, "nuc") - 1)/2.0

def log_neg(ab):
    aTb = partial_transpose(ab,0)
    trace_norm = LA.norm(aTb, "nuc")
    if trace_norm <= 0:
        trace_norm = 1
    return np.log2(trace_norm)

def isEntangled(ab):
    walue, vec = LA.eig(partial_transpose(ab,0))
    if np.sum(np.isreal(walue)) != np.sum(np.shape(walue)):
        return 1
    
    for i in walue:
        if i < 0:
            return 1

    return 0

# print(isEntangled(rho_ab))



def rho_psi(theta):
    psi = np.array([[np.cos(theta/2),0,0, np.sin(theta/2)]])
    return np.kron(psi,psi.T)

rho_ab = rho_psi(np.pi/1.9)

p = 0.1
ab = rho_ab
# print(rho_ab)

def damping_channel_a(ab,p):

    m0 = np.array([[1 , 0], [ 0,np.sqrt(1-p)]])
    m1 = np.array([[0,np.sqrt(p)], [0,0]])

    M0 = np.kron(m0,np.eye(2))
    # print(M0)
    M1 = np.kron(m1,np.eye(2))

    return np.add(np.matmul(np.matmul(M0,ab),M0.T) , np.matmul(np.matmul(M1,ab),M1.T))

def damping_channel_b(ab,p):

    m0 = np.array([[1 , 0], [ 0,np.sqrt(1-p)]])
    m1 = np.array([[0,np.sqrt(p)], [0,0]])

    M0 = np.kron(np.eye(2),m0)
    # print(M0)
    M1 = np.kron(np.eye(2),m1)

    return np.add(np.matmul(np.matmul(M0,ab),M0.T) , np.matmul(np.matmul(M1,ab),M1.T))

def double_damping_channel(ab,p):

    m0 = np.array([[1 , 0], [ 0,np.sqrt(1-p)]])
    m1 = np.array([[0,np.sqrt(p)], [0,0]])

    M00 = np.kron(m0,m0)
    M01 = np.kron(m0,m1)
    M10 = np.kron(m1,m0)
    M11 = np.kron(m1,m1)

    n0 = np.add(np.matmul(np.matmul(M00,ab),M00.T) , np.matmul(np.matmul(M01,ab),M01.T))
    n1 = np.add (np.matmul(np.matmul(M10,ab),M10.T), np.matmul(np.matmul(M11,ab),M11.T))
    return np.add(n0,n1)

def double_damping_channel_1(ab,p):

    m0 = np.array([[1 , 0], [ 0,np.sqrt(1-p)]])
    m1 = np.array([[0,np.sqrt(p)], [0,0]])

    M0 = np.kron(m0,m0)
    # M01 = np.kron(m0,m1)
    # M10 = np.kron(m1,m0)
    M1 = np.kron(m1,m1)

    return np.add(np.matmul(np.matmul(M0,ab),M0.T) , np.matmul(np.matmul(M1,ab),M1.T))


# LN = []
# T = []
# P = []

def Log_Neg_Z(theta, p):
    rho_ab = rho_psi(theta)
    ln = log_neg(damping_channel_b(damping_channel_a(rho_ab,p),p))
    return ln

def Log_Neg_single_damping(theta, p):
    rho_ab = rho_psi(theta)
    ln = log_neg(damping_channel_a(rho_ab,p))
    return ln

t_d = np.linspace(0,np.pi,num=100,endpoint=False)
p_d = np.linspace(0,1,num=100,endpoint=False)

T, P = np.meshgrid(t_d,p_d)
# Z = Log_Neg_Z(T,P)
# Z = np.eye(100)

# for i in range(100):
#     for j in range(100):
#         Z[i,j] = Log_Neg_single_damping(T[i,j],P[i,j])

# for w in range(0,31):
#     for p_i in range(0,10):
#         t = w/10.0
#         p = p_i/10.0
#         rho_ab = rho_psi(t)
#         ln = log_neg(damping_channel_b(damping_channel_a(rho_ab,p),p))
#         # ln = log_neg(double_damping_channel_1(rho_ab,p))
#         LN.append(ln)
#         P.append(p)
#         T.append(t)

# plt.plot(LN)
# plt.show()

# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=plt.figaspect(0.5))


#===============
#  First subplot
#===============
# set up the axes for the first plot
ax = fig.add_subplot(1, 2, 1, projection='3d')

# plot a 3D surface like in the example mplot3d/surface3d_demo
X = T
Y = P

Z = np.eye(100)

for i in range(100):
    for j in range(100):
        Z[i,j] = Log_Neg_single_damping(X[i,j],Y[i,j])

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlabel('θ')
ax.set_ylabel('p')
ax.set_zlabel('Logarithmic Negativity')
ax.set_title('One Channel Noise')

#===============
# Second subplot
#===============
# set up the axes for the second plot
ax = fig.add_subplot(1, 2, 2, projection='3d')

X = T
Y = P
Z = np.eye(100)

for i in range(100):
    for j in range(100):
        Z[i,j] = Log_Neg_Z(X[i,j],Y[i,j])

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.summer_r,
                       linewidth=0, antialiased=False)

ax.set_xlabel('θ')
ax.set_ylabel('p')
ax.set_zlabel('Logarithmic Negativity')
ax.set_title('Two Channel Noise')


plt.show()


# ax = fig.add_subplot(111, projection='3d')


# ax.set_xlabel('θ')
# ax.set_ylabel('p')
# ax.set_zlabel('Logarithmic Negativity')
# Grab some test data.
# X, Y, Z = T,P,LN #axes3d.get_test_data(0.05)

# Plot a basic wireframe.
# ax.plot_wireframe(T, P, Z, rstride=10, cstride=10)

# X = T
# Y = P
# ax.scatter3D(X, Y, Z)
# ax.contour3D(X, Y, Z, 50, cmap='coolwarm')
# ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='nipy_spectral',edgecolor='none')


plt.show()











