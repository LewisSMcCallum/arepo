#!/bin/bash
#SBATCH -J arepo
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=output
#SBATCH --error=error




source /opt/rh/devtoolset-8/enable
export LD_LIBRARY_PATH=/gpfs1/apps/software/devts8/lib:$LD_LIBRARY_PATH
export PATH=/gpfs1/apps/software/devts8/mvapich2/bin:$PATH
export LD_LIBRARY_PATH=/gpfs1/apps/software/devts8/mvapich2/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/gpfs1/apps/conda/lm261/conda/lib:$LD_LIBRARY_PATH
ulimit -s unlimited


mpirun -n $SLURM_NPROCS ../Arepo param.txt 










