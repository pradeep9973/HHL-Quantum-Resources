# author : pradeep 
# hhl disorder and success probability computation
# uniform error

from mpi4py import MPI
import HHL_negativity
# import GGM_module
import time
import random
from itertools import combinations
from qutip.states import *
from qutip.tensor import *
from qutip.qobj import Qobj, isket
import qutip
import numpy as np
from numpy import linalg as LA
import disorder
import sys

def disorder_GGM_PSI2_data_files(range_uniform):
    # start_time = time.time()
    
    file = open("uniform_psi2_test_new" + str(range_uniform)[-1] + ".txt", "w")
    #====================STEPS=====================#
    steps = 10
    b_0_sq = [(i)/(steps*1.0) for i in range(1,steps)]
    b_0 = [i**(0.5) for i in b_0_sq]
    # sigma iteration
    s = range_uniform
    # print(b_0)

    for b in b_0:
        
        # ggm_0 = HHL_negativity.New_GGM(psi_0)
        # ggm_og.append(ggm_0)
        
        #psi3_og

        psi3_og = disorder.PSI_3(b)
        
        #===================MPI==================================
        # MPI Initialization
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Define the number of runs and calculate the number of runs per process
        runs = 5
        runs_per_process = runs // size

        # Initialize the variables for accumulation
        ggm_sum = 0.0
        distance_sum = 0.0
        success_prob_sum = 0.0

        # Set the random seed based on the rank to ensure different random numbers for each process
        random.seed(rank + b*10000)
        print("rank : ", rank)

        # Perform the distributed loop
        for j in range(runs_per_process):
            rn_1 = random.uniform(-s,s) #random.gauss(0, s)
            # print("rn_1 is : ", rn_1, " and rank is  :", rank)
            rn_2 = random.uniform(-s,s) #random.gauss(0, s)
            # print("difference is : ", rn_2 - rn_1, " and rank is  :", rank)
            psi2_d = HHL_negativity.PSI_2_GGM_disorder(b, w=rn_1, u=rn_2)
            psi3_disorder = disorder.PSI_3_disorder(b, w=rn_1, u=rn_2)
            
            ggm_sum += HHL_negativity.New_GGM(psi2_d)/runs/1.0
            # print("ggm_sum : ", ggm_sum, "rank is : ", rank )
            distance_sum += disorder.succ_distance(psi3_og, psi3_disorder)/runs/1.0
            success_prob_sum += disorder.success_prob(psi3_disorder)/runs/1.0

        # Perform reduction operations to gather the accumulated results from all processes
        ggm_average = comm.reduce(ggm_sum, op=MPI.SUM, root=0)
        if rank == 0:
            print("ggm : ", ggm_average, "rank is : ", rank)
        distance_average = comm.reduce(distance_sum, op=MPI.SUM, root=0)
        success_prob_average = comm.reduce(success_prob_sum, op=MPI.SUM, root=0)

        # ggm_sum_total = comm.reduce(ggm_sum, op=MPI.SUM, root=0)
        # distance_sum_total = comm.reduce(distance_sum, op=MPI.SUM, root=0)
        # success_prob_sum_total = comm.reduce(success_prob_sum, op=MPI.SUM, root=0)

        # if rank == 0:
        #     ggm_average = ggm_sum_total / runs
        #     distance_average = distance_sum_total / runs
        #     success_prob_average = success_prob_sum_total / runs

        #=====================MPI-code-over=======================
        # ggm_sum = 0
        # distance_sum = 0
        # success_prob_sum = 0
        # #=================RUNS======================#
        # runs = 100
        # for j in range(runs):
        #     rn_1 = random.gauss(0, s)
        #     rn_2 = random.gauss(0, s)
        #     psi2_d = HHL_negativity.PSI_2_GGM_disorder(b,w = rn_1, u = rn_2)
        #     psi3_disorder = disorder.PSI_3_disorder(b,w = rn_1, u =rn_2)
            
        #     ggm_sum+= HHL_negativity.New_GGM(psi2_d)
        #     distance_sum+= disorder.succ_distance(psi3_og,psi3_disorder)
        #     success_prob_sum+= disorder.success_prob(psi3_disorder)

        # Calculate averages only on the root process
        

        if rank == 0:
            data_write = str(b**2) + '\t' + str(ggm_average) + '\t' + str(distance_average) + '\t' + str(success_prob_average) + '\t' + '\n'
            # print(data_write)
            file.write(data_write)
            # print(b, " done")

    file.close()
    # end_time = time.time()
    # print(end_time-start_time)
    return

if __name__ == "__main__":
    uni_width = float(sys.argv[1])
    print(uni_width)
    disorder_GGM_PSI2_data_files(uni_width)
