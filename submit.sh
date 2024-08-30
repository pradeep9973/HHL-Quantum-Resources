submit.sh 

#!/bin/bash
#PBS 
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -q workq
#PBS -joe
#PBS -V
#PBS -o out.o
#PBS -e out.e


module load intel/2018
module load openmpi/4.0.0
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > pbs_nodes
echo Working directory is $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`
mpirun -np $NPROCS --machinefile $PBS_NODEFILE /home/pradeep/anaconda3/bin/python3 QIC_codes/hhl_disorder.py


