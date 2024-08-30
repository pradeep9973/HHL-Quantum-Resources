# plots

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
import hhl_entanglement


def get_data(file_name):
    #returns a list

    with open(file_name,"r") as file:
        lines = file.readlines()

    data = []
    for i in lines:
        data.append(float(i[:-1]))
    
    return data


# kappa plots


def plot_complexity2d():
    #b_0 =0.3
    
    #get data
    kappa = get_data("Data_files/kappa/kappa.txt")
    coh_r_kappa_norm = [i*2 for i in get_data("Data_files/kappa/coh_r_kappa_norm.txt")]  #while calculating norm i divided by 2.0
    coh_full_kappa_norm = get_data("Data_files/kappa/coh_full_kappa.txt")
    neg_kappa = get_data("Data_files/kappa/neg_kappa.txt")
    ggm_kappa = get_data("Data_files/kappa/ggm_kappa.txt")

    coh_full_norm = []
    for i in range(len(kappa)):
        coh_full_norm.append(coh_full_kappa_norm[i]/(2**(int(np.log2(kappa[i]))+1.0) *4 ))

    #plots
    plt.rcParams["text.usetex"] = True
    fig,ax = plt.subplots()
    ax.plot(kappa, coh_r_kappa_norm, label = r"$\mathcal{C}_{\mathbf{R}}$",linestyle="--",linewidth=2.5)
    ax.plot(kappa, coh_full_norm, label = r"$\mathcal{C}_{\mathbf{\Lambda U R}}$",linestyle="-.",linewidth=2.5)
    ax.plot(kappa, neg_kappa, label = r"$\mathcal{LN}_{\mathbf{\Lambda U}}$",linestyle=":",linewidth=3.5)
    ax.plot(kappa, ggm_kappa, label = r"$\mathcal{E}$",linestyle="-",linewidth=2.5)   
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_xlabel(r"$\bf{\kappa}$",fontsize = 18)
    # ax.set_ylabel(r"\bf{Entanglement measure}",fontsize = 18)
    x_ticks = plt.gca().get_xticks()
    y_ticks = plt.gca().get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # Set the modified tick labels back to the plot
    plt.gca().set_xticklabels(x_tick_labels, fontsize=18)
    plt.gca().set_yticklabels(y_tick_labels, fontsize=18)

    plt.legend(fontsize = 18)
    plt.tight_layout()
    plt.grid()
    plt.show()

    return



def plot_complexity3d():
    
    #b_0 = 0.2, b_1 = 0.3

    #get data
    kappa = get_data("Data_files/kappa/3d_kappa/3d_kappa.txt")
    coh_r_kappa = get_data("Data_files/kappa/3d_kappa/Cr3d_kappa.txt")  #while calculating norm i divided by 2.0
    coh_full_kappa = get_data("Data_files/kappa/3d_kappa/Cfull3d_kappa.txt")
    neg_kappa = get_data("Data_files/kappa/3d_kappa/lneg3d_kappa.txt")
    ggm_kappa = get_data("Data_files/kappa/3d_kappa/ggm3d_kappa.txt")

    coh_full_norm = []
    for i in range(len(kappa)):
        coh_full_norm.append(coh_full_kappa[i]/(2**(int(np.log2(kappa[i]))+1.0) * 6 ))


    #plots
    plt.rcParams["text.usetex"] = True
    fig,ax = plt.subplots()
    ax.plot(kappa, coh_r_kappa, label = r"$\mathcal{C}_R$",linestyle="--")
    ax.plot(kappa, coh_full_norm, label = r"$\mathcal{C}_{\Lambda U R}$",linestyle="-.")
    ax.plot(kappa, neg_kappa, label = r"$\mathcal{LN}_{\Lambda U}$",linestyle=":")
    ax.plot(kappa, ggm_kappa, label = r"$\mathcal{E}$",linestyle="-")   
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_xlabel(r"$\bf{\kappa}$",fontsize = 18)
    # ax.set_ylabel(r"\bf{Entanglement measure}",fontsize = 18)
    x_ticks = plt.gca().get_xticks()
    y_ticks = plt.gca().get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # Set the modified tick labels back to the plot
    plt.gca().set_xticklabels(x_tick_labels, fontsize=18)
    plt.gca().set_yticklabels(y_tick_labels, fontsize=18)

    plt.legend(fontsize = 18)
    plt.tight_layout()
    plt.grid(True)
    plt.show()

    return


def plot_complexity2d3d():
    #b_0 =0.3
    #get data
    kappa = get_data("Data_files/kappa/kappa.txt")
    coh_r_kappa_norm = [i*2 for i in get_data("Data_files/kappa/coh_r_kappa_norm.txt")]  #while calculating norm i divided by 2.0
    coh_full_kappa_norm = get_data("Data_files/kappa/coh_full_kappa.txt")
    neg_kappa = get_data("Data_files/kappa/neg_kappa.txt")
    ggm_kappa = get_data("Data_files/kappa/ggm_kappa.txt")

    coh_full_norm = []
    for i in range(len(kappa)):
        coh_full_norm.append(coh_full_kappa_norm[i]/(2**(int(np.log2(kappa[i]))+1.0) *4 ))

    #plots
    plt.rcParams["text.usetex"] = True
    fig,ax = plt.subplots()
    ax.plot(kappa, coh_r_kappa_norm, label = r"$\mathcal{C}_R^{2d}$",linestyle="--")
    ax.plot(kappa, coh_full_norm, label = r"$\mathcal{C}_{\Lambda U R}^{2d}$",linestyle="-.")
    ax.plot(kappa, neg_kappa, label = r"$\mathcal{LN}_{\Lambda U}^{2d}$",linestyle=":")
    ax.plot(kappa, ggm_kappa, label = r"$\mathcal{E}^{2d}$",linestyle="-")   
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_xlabel(r"$\bf{\kappa}$",fontsize = 18)
    # ax.set_ylabel(r"\bf{Entanglement measure}",fontsize = 18)
    # x_ticks = plt.gca().get_xticks()
    # y_ticks = plt.gca().get_yticks()
    # x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    # y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # # Set the modified tick labels back to the plot
    # plt.gca().set_xticklabels(x_tick_labels, fontsize=18)
    # plt.gca().set_yticklabels(y_tick_labels, fontsize=18)

    # plt.legend(fontsize = 18)
    # plt.tight_layout()
    # plt.grid()
    # plt.show()
    #b_0 = 0.2, b_1 = 0.3

    #get data
    kappa = get_data("Data_files/kappa/3d_kappa/3d_kappa.txt")
    coh_r_kappa = get_data("Data_files/kappa/3d_kappa/Cr3d_kappa.txt")  #while calculating norm i divided by 2.0
    coh_full_kappa = get_data("Data_files/kappa/3d_kappa/Cfull3d_kappa.txt")
    neg_kappa = get_data("Data_files/kappa/3d_kappa/lneg3d_kappa.txt")
    ggm_kappa = get_data("Data_files/kappa/3d_kappa/ggm3d_kappa.txt")

    coh_full_norm = []
    for i in range(len(kappa)):
        coh_full_norm.append(coh_full_kappa[i]/(2**(int(np.log2(kappa[i]))+1.0) * 6 ))


    #plots
    # plt.rcParams["text.usetex"] = True
    # fig,ax = plt.subplots()
    ax.plot(kappa, coh_r_kappa, label = r"$\mathcal{C}_R^{3d}$",linestyle="-",marker="x",markersize =3)
    ax.plot(kappa, coh_full_norm, label = r"$\mathcal{C}_{\Lambda U R}^{3d}$",linestyle=":")
    ax.plot(kappa, neg_kappa, label = r"$\mathcal{LN}_{\Lambda U}^{3d}$",linestyle="-",marker="d",markersize =3)
    ax.plot(kappa, ggm_kappa, label = r"$\mathcal{E}^{3d}$",linestyle="--")   
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    # ax.set_xlabel(r"$\bf{\kappa}$",fontsize = 18)
    # ax.set_ylabel(r"\bf{Entanglement measure}",fontsize = 18)
    x_ticks = plt.gca().get_xticks()
    y_ticks = plt.gca().get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # Set the modified tick labels back to the plot
    plt.gca().set_xticklabels(x_tick_labels, fontsize=18)
    plt.gca().set_yticklabels(y_tick_labels, fontsize=18)

    plt.legend(fontsize = 11)
    plt.tight_layout()
    plt.grid(True)
    plt.show()

    return



# 2D plots
def plot_2d_ent():
       
    #get data
    b0_sq = get_data("Data_files/2D/b0_sq.txt")
    ggm = get_data("Data_files/2D/ggm2d.txt")
    lneg = get_data("Data_files/2D/lneg2d.txt")
    C_r = get_data("Data_files/2D/C_r_2d.txt")
    C_full = get_data("Data_files/2D/C_full_2d.txt")
    sp = get_data("Data_files/2D/sp.txt")


    #plots
    plt.rcParams["text.usetex"] = True
    fig,ax = plt.subplots()
    ax.plot(b0_sq, ggm, label = r"$\mathcal{E}$",linestyle="-",linewidth=2.5)
    ax.plot(b0_sq, lneg, label = r"$\mathcal{LN}_{\mathbf{\Lambda U}}$",linestyle=":",linewidth=3.5)
    ax.plot(b0_sq, C_r, label = r"$\mathcal{C}_{\mathbf{R}}$",linestyle="-.",linewidth=2.5)
    ax.plot(b0_sq, C_full, label = r"$\mathcal{C}_{\mathbf{\Lambda U R}}$",linestyle="--",linewidth=2.5) 
    ax.plot(b0_sq,sp, label = r"$\mathcal{SP}$", linestyle ="-",  marker = 'o',markersize = 3,linewidth=2.5)  
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.set_xlabel(r"$\bf{b_0^2}$",fontsize = 20,fontweight="bold")
    # x_tick_labels = [r"$\bf{" + tick.get_text() + r"}$" for tick in plt.gca().get_xticklabels()]
    # ax.set_xticklabels(x_tick_labels,fontdict = {"fontweight":"bold"})


    # Get the tick labels for the x and y axes
    # Get the current tick locations and labels for both x and y axes
    x_ticks = plt.gca().get_xticks()
    y_ticks = plt.gca().get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # Set the modified tick labels back to the plot
    plt.gca().set_xticklabels(x_tick_labels, fontsize=18)
    plt.gca().set_yticklabels(y_tick_labels, fontsize=18)



    # ax.set_ylabel(r"\bf{Entanglement measure}",fontsize = 20, rotation=90)
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return



# 3d plots

def read3d_data(file_name):
    ggm_3d = np.zeros((50, 50))
    with open(file_name, "r") as file:
        lines = file.readlines()
        for line in lines:
            values = line.strip().split()
            i = int(float(values[0]))
            j = int(float(values[1]))
            ggm_ij = float(values[2])
            ggm_3d[i, j] = ggm_ij
    return ggm_3d


def plot_3d_data(data_3d): 

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
    contourf_plot = plt.contourf(X, Y, Z,  levels=10,cmap='RdBu')
    contour_labels = plt.contour(X, Y, Z, levels=10, colors='black', linewidths=0.7)
    # plt.clabel(contour_labels, inline=True, fontsize=8)

    # Set labels and title
    plt.xlabel(r"$\mathbf{b_{0}^{2}}$", fontsize=18)
    plt.ylabel(r"$\mathbf{b_{1}^{2}}$", fontsize=18, rotation=90)
    plt.xticks(fontsize=18)  # Change the font size of x-axis tick labels
    plt.yticks(fontsize=18)  # Change the font size of y-axis tick labels
    plt.clabel(contour_labels, inline=True, fontsize=18)
    # plt.text(1.05, 1.05, r'$\mathcal{C}_{\Lambda U R}$', fontsize=16, transform=plt.gcf().transFigure)

    x_ticks = plt.gca().get_xticks()
    y_ticks = plt.gca().get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # Set the modified tick labels back to the plot
    plt.gca().set_xticklabels(x_tick_labels, fontsize=20)
    plt.gca().set_yticklabels(y_tick_labels, fontsize=20)


    #title = r"$\mathcal{SP}$"
    # plt.title(title)
    colorbar = plt.colorbar(contourf_plot)
    # Increase the tick size of the colorbar
    colorbar.ax.tick_params(labelsize=18)
    # colorbar.set_label(r'$\mathcal{C}_{\Lambda U R}$', fontsize=14, labelpad=15)
    # Create a new axis on top of the colorbar axis to display the label
    colorbar_label = r'$\mathcal{SP}$'

    # Specify the desired colorbar tick positions
    num_ticks = len(colorbar.get_ticks())
    num_colorbar_ticks = num_ticks
    colorbar_ticks = np.linspace(colorbar.vmin, colorbar.vmax, num_colorbar_ticks)
    
    # Get the colorbar tick labels and modify them to be bold using LaTeX
    colorbar_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in colorbar_ticks]

    # Set the modified colorbar tick positions and labels
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels(colorbar_tick_labels, fontsize=20, weight='bold')




    colorbar_coords = colorbar.ax.get_position()
    plt.annotate(colorbar_label, xy=(colorbar_coords.x0 + 0.7 * colorbar_coords.width, colorbar_coords.y1 + 0.04), xycoords='figure fraction', ha='center', fontsize=16)
    plt.tight_layout()


    # Show the plot
    plt.savefig("pdf_images/sp3d.pdf")
    # plt.show()
    return


def plot_file(file_name):
    data = read3d_data(file_name)
    plot_3d_data(data)
    return


def plot_3d_data_2plots(data_3d, ax, ano_x,ano_y, annotation ,cmap_color):
    plt.rcParams["text.usetex"] = True
    x = np.arange(0, data_3d.shape[1])
    y = np.arange(0, data_3d.shape[0])
    x_sq = (x/100.0)
    y_sq = (y/100.0)
    X, Y = np.meshgrid(x_sq, y_sq)
    Z = data_3d

    contourf_plot = ax.contourf(X, Y, Z, levels=10, cmap=cmap_color) 
    contour_labels = ax.contour(X, Y, Z, levels=10, colors='black', linewidths=0.7)
    # ax.clabel(contour_labels, inline=True, fontsize=8)

    ax.set_xlabel(r"$\mathbf{b_{0}^{2}}$", fontsize=22)
    ax.set_ylabel(r"$\mathbf{b_{1}^{2}}$", fontsize=22, rotation=90)
    ax.tick_params(axis='both', labelsize=22)
    ax.clabel(contour_labels, inline=True, fontsize=22)



    x_ticks = ax.get_xticks()
    y_ticks = ax.get_yticks()
    x_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in x_ticks]
    y_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in y_ticks]

    # Set the modified tick labels back to the plot
    ax.set_xticklabels(x_tick_labels, fontsize=22)
    ax.set_yticklabels(y_tick_labels, fontsize=22)

    colorbar = plt.colorbar(contourf_plot)
    colorbar.ax.tick_params(labelsize=22)
    # Get the colorbar tick labels and modify them to be bold using LaTeX

    # Specify the desired colorbar tick positions
    num_ticks = len(colorbar.get_ticks())
    num_colorbar_ticks = num_ticks
    colorbar_ticks = np.linspace(colorbar.vmin, colorbar.vmax, num_colorbar_ticks)
    
    # Get the colorbar tick labels and modify them to be bold using LaTeX
    colorbar_tick_labels = [r"$\mathbf{" + f"{tick:.2f}" + r"}$" for tick in colorbar_ticks]

    # Set the modified colorbar tick positions and labels
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels(colorbar_tick_labels, fontsize=22, weight='bold')


    colorbar_label = annotation #r'$\mathcal{C}_{\Lambda U R}$'
    colorbar_coords = colorbar.ax.get_position()
    annotation_x = colorbar_coords.x0 + ano_x * colorbar_coords.width
    annotation_y = colorbar_coords.y1 + ano_y  # Adjusting the vertical alignment

    ax.annotate(colorbar_label, xy=(annotation_x, annotation_y),
                xycoords='figure fraction', ha='center', va='center', fontsize=22)



def plot_3d_data_side_by_side(data_1, data_2, cmap_color = "cool"):
    plt.rcParams["text.usetex"] = True
    fig, axes = plt.subplots(2, 1, figsize=(8, 12))
    # plot_3d_data_2plots(data_1, axes[0], 0.7, 0.06, r'$\mathcal{C}_{\Lambda U R}$',cmap_color)
    # plot_3d_data_2plots(data_2, axes[1], 0.7, 0.012, r'$\mathcal{C}_{R}$',cmap_color)
    plot_3d_data_2plots(data_1, axes[0], 0.7, 0.06, r'$\mathcal{E}$',cmap_color)
    plot_3d_data_2plots(data_2, axes[1], 0.7, 0.012, r'$\mathcal{LN}_{|Lambda U}$',cmap_color)
    plt.tight_layout()
    plt.savefig("pdf_images/ggm_lneg_3d_b.pdf")
    # plt.show()

def plot_3d_ggm_lneg():#data_1, data_2, cmap_color = "cool"):
    data_1 = read3d_data("Data_files/3D/ggm3d.txt")
    # data_2 = read3d_data("Data_files/3D/lneg3d.txt")
    data_2 = read3d_data("Data_files/3D/lneg.txt")

    cmap_color  = "RdBu"
    plt.rcParams["text.usetex"] = True
    fig, axes = plt.subplots(2, 1, figsize=(8, 12))
    # plot_3d_data_2plots(data_1, axes[0], 0.7, 0.06, r'$\mathcal{C}_{\Lambda U R}$',cmap_color)
    # plot_3d_data_2plots(data_2, axes[1], 0.7, 0.012, r'$\mathcal{C}_{R}$',cmap_color)
    plot_3d_data_2plots(data_1, axes[0], 0.47, 0.06, r'$a) \;\mathbf{\mathcal{E}}$',cmap_color)
    plot_3d_data_2plots(data_2, axes[1], 0.55, 0.012, r'$b) \;\mathbf{\mathcal{LN}_{\mathbf{\Lambda U}}}$',cmap_color)
    plt.tight_layout()
    plt.savefig("pdf_images/ggm_lneg_3d_ln2.pdf")
    # plt.show()

def plot_3d_coherence():#data_1, data_2, cmap_color = "cool"):
    data_1 = read3d_data("Data_files/3D/C_full_norm.txt")
    data_2 = read3d_data("Data_files/3D/C_r_norm.txt")
    cmap_color  = "RdBu"
    plt.rcParams["text.usetex"] = True
    fig, axes = plt.subplots(2, 1, figsize=(8, 12))
    plot_3d_data_2plots(data_1, axes[0], 0.55, 0.06, r'$a) \; \mathbf{\mathcal{C}_{\mathbf{\Lambda U R}}}$',cmap_color)
    plot_3d_data_2plots(data_2, axes[1], 0.55, 0.012, r'$b) \; \mathbf{\mathcal{C}_{\mathbf{R}}}$',cmap_color)
    # plot_3d_data_2plots(data_1, axes[0], 0.7, 0.06, r'$\mathcal{E}$',cmap_color)
    # plot_3d_data_2plots(data_2, axes[1], 0.7, 0.012, r'$\mathcal{LN}_{\Lambda U}$',cmap_color)
    # axes[0].set_title("a)")
    plt.tight_layout()
    # plt.show()

    plt.savefig("pdf_images/coherence_3d.pdf")
