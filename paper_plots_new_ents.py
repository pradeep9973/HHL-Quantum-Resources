#author : pradeep kumar
# plots

import GGM_module
import time
import random
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
import disorder
import csv

def ggm_plot():
    b0 = [i/1000.0 for i in range(1001)]
    b0_sq = [i**2 for i in b0]

    psi_list = [HHL_negativity.PSI_GGM(i, 1, 2) for i in b0]
    psi_3_list = [HHL_negativity.PSI_3(i, 1, 2) for i in b0]
    ggm = [HHL_negativity.New_GGM(i) for i in psi_list]


    efficiency = [disorder.success_prob(i) for i in psi_3_list]


    # psi_list_1 = [HHL_negativity.PSI_GGM(i, 2, 5) for i in b0]
    # ggm_1 = [HHL_negativity.New_GGM(i) for i in psi_list]

    plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots(figsize=(6, 3.6))
    x = b0_sq
    y1 = ggm
    y2 = efficiency
    line1 = ax.plot(x, y1, label=r'$\mathcal{E}$',color = "black")
    ax2 = ax.twinx()
    line2 = ax2.plot(x,y2, label=r"Success Probability", color = "#26718e", linestyle='dashed')
    
    # ax.plot(x, y2, label=r'lambda_{23}')

    # Add grid and legend
    # ax.grid(True)
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    plt.legend(lines, labels)
    # ax.legend()
    # ax2.legend()

    # Set title and axes labels in LaTeX
    # ax.set_title(r'\textbf{Entanglement in HHL Algorithm}')
    ax.set_xlabel(r"$b_{0}^{2}$",  fontsize=14)
    ax.set_ylabel(r'$\mathcal{E}$',  fontsize=14, color="black")
    ax2.set_ylabel(r"Success Probability", fontsize=14, color = "#26718e")
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12, colors="black")
    ax2.tick_params(axis='y', labelsize=12, colors="#26718e")
    fig.tight_layout()

    # Save plot in EPS format
    fig.savefig("fig_GGM_SP.eps", format="eps")
    # plt.show()

def negativity_plot():
    ratio = [i / 100.0 for i in range(2, 100)]
    psi_list = [HHL_negativity.PSI_2(0.1, 100, i) for i in range(2, 100)]
    counter = 0
    negativity = []
    for i in psi_list:
         qubits = int(np.log2(i.shape[0]))
         sys = [i for i in range(qubits)]
         sub_sys = [0 for i in range(qubits)]
         sub_sys[-2] = 1
         neg = HHL_negativity.custom_negativity(i, sys[:-1], sub_sys[:-1])
         negativity.append(neg)
         print(counter)
         counter += 1
    
    plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots()
    x = ratio
    y1 = negativity
    # y2 = ggm_1
    ax.plot(x, y1, label=r'Negativity',color='#26718e')
    # ax.plot(x, y2, label=r'lambda_{23}')

    # Add grid and legend
    # ax.grid(True)
    ax.legend()

    # Set title and axes labels in LaTeX
    ax.set_title(r'\textbf{Entanglement vs ratio of eigenvalues}')
    ax.set_xlabel(r"$\mathrm{\frac{\lambda_{2}}{\lambda_{1}}}$")
    ax.set_ylabel(r'Negativity')

    # Save plot in EPS format
    fig.savefig("fig_2.eps", format="eps")
    # plt.show()

def efficiency_vs_ratio_plot():
    def success_prob(PSI):
        prob = abs(PSI[1])**2 + abs(PSI[3])**2 
        return prob[0][0]
    b = 0.1
    # ratio = [i/10000 for i in range(1,1000)]
    lambda_1 = 40
    lambda_2 = 200
    lambdas = [i for i in range(1,lambda_1)]
    ratio = [i/lambda_2 for i in lambdas]
    actual_prob = []
    for l in lambdas :
        d = []
        actual = success_prob(HHL_negativity.PSI_3(b,l,lambda_2))
        actual_prob.append(actual)
        print(b,"-->", actual)
    psi_list = [HHL_negativity.PSI_GGM(0.2,i,lambda_2) for i in lambdas]
    def print_ggm(psi):
        ggm = HHL_negativity.New_GGM(psi)
        print("ggm == ",ggm)
        return ggm

    ggm =[]
    counter = 0
    for i in psi_list:
        g = HHL_negativity.New_GGM(i)
        ggm.append(g)
        print(counter,g)
        counter+=1
    plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots(figsize=(6, 3.6))
    x = ratio
    y1 = actual_prob
    y2 = ggm
    ax.plot(x, y1, label=r'Efficiency',color='black')
    ax2 = ax.twinx()
    ax2.plot(x, y2, label=r'GGM', color="black", linestyle='dashed')

    # Add grid and legend
    # ax.grid(True)
    ax.legend()

    # Set title and axes labels in LaTeX
    # ax.set_title(r'\textbf{Entanglement vs ratio of eigenvalues}')
    ax.set_xlabel(r"$\mathrm{\frac{\lambda_{2}}{\lambda_{1}}}$", fontsize=13)
    # ax.set_xlabel(r"$b_{0}^{2}$")

    ax.set_ylabel(r'Efficiency',fontsize=12)
    ax2.set_ylabel(r'GGM',fontsize=12),
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax2.tick_params(axis='y', labelsize=12)
    fig.tight_layout()

    # fig.savefig("fig_4.eps", format="eps")
    plt.show()


def disorder_GGM():
    start_time = time.time()
    b_0_sq = [(i)/1000.0 for i in range(1,1000)]
    b_0 = [i**(0.5) for i in b_0_sq]

    # without disorder
    ggm_og =[]
    for b in  b_0:
        psi_0 = disorder.PSI_3_disorder(b,w = 0, u =0)
        ggm_0 = HHL_negativity.New_GGM(psi_0)
        ggm_og.append(ggm_0)

    #disorder
    sigma = [0.1, 0.2, 0.4, 0.8] 

    #plot figure object
    plt.rcParams["text.usetex"] = True
    # fig, ax = plt.subplots(figsize=(6, 3.6))
    fig, ax = plt.subplots(figsize=(10, 6))


    # sigma iteration
    for s in sigma:
        print("+++++++++++++sigma = ",s," +++++++++++++++++") 
        ggm_av =[]  
        for b in b_0:
            # psi_0 = disorder.PSI_3_disorder(b,w = 0, u =0)
            # ggm_0 = HHL_negativity.New_GGM(psi_0)
            # ggm_og.append(ggm_0)
            ggm_sum = 0
            runs = 2
            for j in range(runs):
                rn_1 = random.gauss(0, s)
                rn_2 = random.gauss(0, s)
                psi_d = disorder.PSI_3_disorder(b,w = rn_1, u = rn_2)
                ggm_sum+= HHL_negativity.New_GGM(psi_d)

            ggm_average = ggm_sum/runs
            ggm_av.append(ggm_average)
            # print(int(b*100),ggm_0, ggm_average)
        
        x = b_0_sq
        y = ggm_av
        ax.plot(x, y, label=r'$\sigma = $' + str(s))
        # plt.plot(b_0_sq,ggm_av,label=s)
    
    x = b_0_sq
    y = ggm_og
    ax.plot(x, y, label=r'$\sigma = 0$')

    # Add grid and legend
    ax.grid(True)
    ax.legend()

    # Set title and axes labels in LaTeX
    ax.set_xlabel(r"$b_{0}^{2}$",  fontsize=12)
    ax.set_ylabel(r'GGM',  fontsize=12)
    ax.tick_params(axis='x', labelsize=12)
    # fig.tight_layout()
    end_time = time.time()
    print(end_time-start_time)
    plt.show()
    # plt.savefig("disorder.eps", format="eps")
    return


def disorder_GGM_data_files(sigma):
    # start_time = time.time()
    
    file = open("sigma_" + str(sigma)[-1] + ".txt", "w")

    b_0_sq = [(i)/1000.0 for i in range(1,1000)]
    b_0 = [i**(0.5) for i in b_0_sq]



    # sigma iteration
    s = sigma

    for b in b_0:
        # psi_0 = disorder.PSI_3_disorder(b,w = 0, u =0)
        # ggm_0 = HHL_negativity.New_GGM(psi_0)
        # ggm_og.append(ggm_0)
        ggm_sum = 0
        runs = 10000
        for j in range(runs):
            rn_1 = random.gauss(0, s)
            rn_2 = random.gauss(0, s)
            psi_d = disorder.PSI_3_disorder(b,w = rn_1, u = rn_2)
            ggm_sum+= HHL_negativity.New_GGM(psi_d)

        ggm_average = ggm_sum/runs
        file.write(str(b**2) + '\t' + str(ggm_average) + '\n')


    file.close()
    # end_time = time.time()
    # print(end_time-start_time)
    return


def disorder_GGM_PSI2_data_files(sigma):
    # start_time = time.time()
    
    file = open("sigma_psi2_test" + str(sigma)[-1] + ".txt", "w")

    steps = 20
    b_0_sq = [(i)/(steps*1.0) for i in range(1,steps)]
    b_0 = [i**(0.5) for i in b_0_sq]


    # sigma iteration
    s = sigma

    for b in b_0:
        
        # ggm_0 = HHL_negativity.New_GGM(psi_0)
        # ggm_og.append(ggm_0)
        
        #psi3_og
        psi3_og = disorder.PSI_3(b)
        
        ggm_sum = 0
        distance_sum = 0
        success_prob_sum = 0
        runs = 100
        for j in range(runs):
            rn_1 = random.gauss(0, s)
            rn_2 = random.gauss(0, s)
            psi2_d = HHL_negativity.PSI_2_GGM_disorder(b,w = rn_1, u = rn_2)
            psi3_disorder = disorder.PSI_3_disorder(b,w = rn_1, u =rn_2)
            
            ggm_sum+= HHL_negativity.New_GGM(psi2_d)
            distance_sum+= disorder.succ_distance(psi3_og,psi3_disorder)
            success_prob_sum+= disorder.success_prob(psi3_disorder)

        ggm_average = ggm_sum/runs
        distance_avg = distance_sum/runs
        success_prob_avg = success_prob_sum/runs
        data_write = str(b**2) + '\t' + str(ggm_average) + '\t' + str(distance_avg) + '\t' + str(success_prob_avg) + '\t' + '\n'
        # print(data_write)
        file.write(data_write)
        # print(b, " done")

    file.close()
    # end_time = time.time()
    # print(end_time-start_time)
    return

# 3D plot


# from mpl_toolkits import mplot3d
# from mpl_toolkits.mplot3d import Axes3D
# import random

def disorder_GGM_3D():
    b_0_sq = [(i)/100.0 for i in range(1,100)]
    b_0 = [i**(0.5) for i in b_0_sq]

    # without disorder
    ggm_og =[]
    for b in b_0:
        psi_0 = disorder.PSI_3_disorder(b,w = 0, u =0)
        ggm_0 = HHL_negativity.New_GGM(psi_0)
        ggm_og.append(ggm_0)

    #disorder
    sigma = [0.1, 0.2, 0.4, 0.8] 
    ggm_av = np.zeros((len(b_0), len(sigma)))

    # sigma iteration
    for i, s in enumerate(sigma):
        print("+++++++++++++sigma = ",s," +++++++++++++++++") 
        for j, b in enumerate(b_0):
            ggm_sum = 0
            runs = 2
            for k in range(runs):
                rn_1 = random.gauss(0, s)
                rn_2 = random.gauss(0, s)
                psi_d = disorder.PSI_3_disorder(b,w = rn_1, u = rn_2)
                ggm_sum += HHL_negativity.New_GGM(psi_d)

            ggm_average = ggm_sum/runs
            ggm_av[j, i] = ggm_average
            print(int(b*100), ggm_0, ggm_average)

    # create meshgrid
    B0, SIGMA = np.meshgrid(b_0_sq, sigma)

    # create figure and 3D axes
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # plot surface
    ax.plot_surface(B0, SIGMA, ggm_av, cmap='viridis')

    # Set title and axes labels in LaTeX
    ax.set_xlabel(r"$b_{0}^{2}$", fontsize=12)
    ax.set_ylabel(r'$\sigma$', fontsize=12)
    ax.set_zlabel(r'GGM', fontsize=12)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='z', labelsize=12)

    plt.show()
    return


def test_ggm_disorder_vs_succ_prob():
    return



def disorder_vs_probability():

    # b_0_sq = [(i)/100.0 for i in range(1,100)]
    # b_0 = [i**(0.5) for i in b_0_sq]
    b_0 = 0.2

    # disorder.success_prob(psi)

    #disorder
    # sigma = [0.1, 0.2, 0.4, 0.8] 
    sigma = [i/10.0 for i in range(0,10)] 

    # without disorder
    ggm =[]
    success_prob_sigma = [] 

    #plot figure object
    # plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots(figsize=(6, 3.6))

    # sigma iteration
    for s in sigma:
        print("+++++++++++++sigma = ",s," +++++++++++++++++") 
        ggm_sum = 0
        success_prob_sum = 0
        runs = 100
        for j in range(runs):
            rn_1 = random.gauss(0, s)
            rn_2 = random.gauss(0, s)
            psi_d = disorder.PSI_3_disorder(b_0,w = rn_1, u = rn_2)
            ggm_sum+= HHL_negativity.New_GGM(psi_d)
            success_prob_sum += disorder.success_prob(psi_d)

        ggm_average = ggm_sum/runs
        success_prob_av = success_prob_sum/runs
        ggm.append(ggm_average)
        success_prob_sigma.append(success_prob_av)
        print(s, ggm_average, success_prob_av)
    
    x = sigma
    y1 = ggm
    y2 = success_prob_sigma
    ax.plot(x, y1)
    ax.plot(x, y2)
    

    # Add grid and legend
    ax.grid(True)
    ax.legend()

    # Set title and axes labels in LaTeX
    # ax.set_xlabel(r"$b_{0}^{2}$",  fontsize=12)
    # ax.set_ylabel(r'GGM',  fontsize=12)
    # ax.tick_params(axis='x', labelsize=12)
    # fig.tight_layout()

    plt.show()
    # plt.savefig("disorder.eps", format="eps")
        

    return


def get_sigmas(file_name):
    file = open(file_name,"r")
    lines = file.readlines()
    b_sq =[]
    ggm = []
    sd = []
    sp = []

    for l in lines:
        a = l.split("\t")
        b_sq.append(float(a[0]))
        ggm.append(float(a[1]))
        sd.append(float(a[2]))
        sp.append(float(a[3]))
    return b_sq,ggm,sd,sp

def plot_sigmas(file_names):
    plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for f in file_names:
        b_sq, ggm,sp,sd = get_sigmas(f)
        x = b_sq
        y = ggm
        ax.plot(x, y, label=r'$\sigma$ = 0.'+f[-5:-4])
        # ax.plot(x,sp,label = r'error of $\sigma$ = 0.'+f[-5:-4])

    
    ax.grid(True)
    ax.legend(fontsize=18)
    # Set title and axes labels in LaTeX
    # ax.set_title(r'\textbf{Entanglement in HHL Algorithm}')
    ax.set_xlabel(r"$b_{0}^{2}$",  fontsize=18)
    ax.set_ylabel(r'$\langle \mathcal{E} \rangle $',  fontsize=18, rotation=90)
    # ax2.set_ylabel(r"Efficiency", fontsize=12)
    # ax.tick_params(axis='x', labelsize=20)
    # ax.tick_params(axis='y', labelsize=20)
    x_ticks = plt.gca().get_xticks()
    y_ticks = plt.gca().get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]
    plt.gca().set_xticklabels(x_tick_labels, fontsize=18)
    plt.gca().set_yticklabels(y_tick_labels, fontsize=18)

    fig.tight_layout()

    # Save plot in EPS format
    # fig.savefig("disorder_27-7-23.eps", format="eps")
    plt.show()

def plot_sigmas_ents():
    plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.grid(True)


    # Set title and axes labels in LaTeX
    # ax.set_title(r'\textbf{Entanglement in HHL Algorithm}')
    ax.set_xlabel(r"$b_{0}^{2}$", fontsize=18, weight='bold')
    # ax.set_ylabel(r'$ \langle \mathcal{C}_{\mathbf{R}} \rangle$', fontsize=18, weight='bold')
    ax.set_ylabel(r"$\langle \mathcal{LN}_{\mathbf{\Lambda U}} \rangle$", fontsize=18, rotation=90)

    ax2 = ax.twinx()
    # ax2.set_ylabel(r"$\langle \mathcal{LN}_{\mathbf{\Lambda U}} \rangle$", fontsize=18, rotation=90)
    ax2.set_ylabel(r'$ \langle \mathcal{C}_{\mathbf{R}} \rangle$', fontsize=18, weight='bold')

    # plt.tick_params(axis='both', which='both', labelsize=18, weight='bold')
#    file_names = ["new_data_test/sigma_psi2_2.txt", "new_data_test/sigma_psi2_4.txt"]
#    for f in file_names:
#        b_sq, ggm,sd,sp = get_sigmas(f)
#        x = b_sq
#        y = ggm
#        y1 = sd
#        # y2 = sp
#        ax.plot(x, y1, label=r'$\sigma$ = 0.'+f[-5:-4] )
#        #ax2.plot(x,y)

    # files_entanglement = ["avg_ent/avg_entanglement_0.1.txt", "avg_ent/avg_entanglement_0.2.txt","avg_ent/avg_entanglement_0.4.txt"]
    # files_entanglement = ["avg_ent/avg_entanglement_0.1.txt", "avg_ent/avg_entanglement_0.2.txt"]
    # files_entanglement = ["Data_files/disorder/avg_ent_new1.txt", "Data_files/disorder/avg_ent_new2.txt", "avg_ent/avg_entanglement_0.4.txt"]
    files_entanglement = ["Data_files/disorder/avg_ent_new1.txt","Data_files/disorder/avg_ent_new2.txt","Data_files/disorder/avg_ent_new4.txt"]


    GGM = []
    log_neg =[]
    C_r = []
    for file in files_entanglement:
        f = open(file,"r")
        csv_reader = csv.reader(f)
        b0 = []
        b0_sq = []
        ggm =[]
        cr =[]
        ln = []
        for row in csv_reader:
            b0.append(float(row[0]))
            b0_sq.append(float(row[1]))
            ggm.append(float(row[2]))
            cr.append(float(row[3]))
            ln.append(float(row[4]))
        GGM.append(ggm)
        log_neg.append(ln)
        C_r.append(cr)
        ax.plot(b0_sq,ln, label=r'$\sigma$ = 0.'+ file[-5:-4])
        ax2.plot(b0_sq,cr,linestyle=":")
        f.close()
        # ax.plot(x,y2)
        # ax.plot(x,sp,label = r'error of $\sigma$ = 0.'+f[-5:-4])
    
        # Keep your x-axis tick labels as they are
    x_ticks = ax.get_xticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    ax.set_xticklabels(x_tick_labels, fontsize=18)

    # Modify y-axis tick labels for the first y-axis (ax)
    y_ticks = ax.get_yticks()
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]
    ax.set_yticklabels(y_tick_labels, fontsize=18)

    # Modify y-axis tick labels for the second y-axis (ax2)
    y2_ticks = ax2.get_yticks()
    y2_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y2_ticks]
    ax2.set_yticklabels(y2_tick_labels, fontsize=18)

    num_ticks = 6  # Set the desired number of ticks
    ax.yaxis.set_major_locator(plt.MaxNLocator(num_ticks))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(num_ticks))


    
    ax.legend(fontsize=18)
    fig.tight_layout()

    # Save plot in EPS format
    # fig.savefig("disorder.eps", format="eps")
    plt.show()



def plot_sigmas_ents2():
    plt.rcParams["text.usetex"] = True
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.grid(True)


    # Set title and axes labels in LaTeX
    # ax.set_title(r'\textbf{Entanglement in HHL Algorithm}')
    ax.set_xlabel(r"$b_{0}^{2}$", fontsize=18, weight='bold')
    # ax.set_ylabel(r'$ \langle \mathcal{C}_R \rangle$', fontsize=18, weight='bold')

    # ax2 = ax.twinx()
    # ax2.set_ylabel(r"$\langle \mathcal{LN} \rangle$", fontsize=18, rotation=90)

    # plt.tick_params(axis='both', which='both', labelsize=18, weight='bold')
#    file_names = ["new_data_test/sigma_psi2_2.txt", "new_data_test/sigma_psi2_4.txt"]
#    for f in file_names:
#        b_sq, ggm,sd,sp = get_sigmas(f)
#        x = b_sq
#        y = ggm
#        y1 = sd
#        # y2 = sp
#        ax.plot(x, y1, label=r'$\sigma$ = 0.'+f[-5:-4] )
#        #ax2.plot(x,y)

    col =["red", "blue", "green"]
    count = 0
    files_entanglement = ["Data_files/disorder/avg_ent_new1.txt", "avg_ent/avg_entanglement_0.1.txt", "avg_ent/avg_entanglement_0.2.txt"]
    GGM = []
    log_neg =[]
    C_r = []
    for file in files_entanglement:
        f = open(file,"r")
        csv_reader = csv.reader(f)
        b0 = []
        b0_sq = []
        ggm =[]
        cr =[]
        ln = []
        for row in csv_reader:
            b0.append(float(row[0]))
            b0_sq.append(float(row[1]))
            ggm.append(float(row[2]))
            cr.append(float(row[3]))
            ln.append(float(row[4]))
        GGM.append(ggm)
        log_neg.append(ln)
        C_r.append(cr)
        ax.plot(b0_sq,cr, label=r'$\sigma$ = 0.'+file[-5:-4],color=col[count])
        ax.plot(b0_sq,ln,color=col[count] )
        count+=1
        f.close()
        # ax.plot(x,y2)
        # ax.plot(x,sp,label = r'error of $\sigma$ = 0.'+f[-5:-4])
    
        # Keep your x-axis tick labels as they are
    x_ticks = ax.get_xticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    ax.set_xticklabels(x_tick_labels, fontsize=18)

    # Modify y-axis tick labels for the first y-axis (ax)
    y_ticks = ax.get_yticks()
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]
    ax.set_yticklabels(y_tick_labels, fontsize=18)

    # # Modify y-axis tick labels for the second y-axis (ax2)
    # y2_ticks = ax2.get_yticks()
    # y2_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y2_ticks]
    # ax2.set_yticklabels(y2_tick_labels, fontsize=18)

    
    ax.legend(fontsize=12)
    fig.tight_layout()

    # Save plot in EPS format
    # fig.savefig("disorder.eps", format="eps")
    plt.show()



if __name__ == "__main__":
    #disorder_GGM_PSI2_data_files(0.4)
    # plot_sigmas(["new_data_tsest/sigma_psi2_0.txt","new_data_test/sigma_psi2_2.txt", "new_data_test/sigma_psi2_4.txt", "new_data_test/sigma_psi2_8.txt"])
    #plot_sigmas_ggm_error(["new_data_test/sigma_psi2_2.txt", "new_data_test/sigma_psi2_4.txt"])
    plot_sigmas_ents()
    # pass