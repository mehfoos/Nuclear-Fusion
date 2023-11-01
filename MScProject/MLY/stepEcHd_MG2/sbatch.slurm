#!/bin/bash
#SBATCH --job-name=stepEcHd_MG2
#SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mly509@york.ac.uk      # Where to send mail
#SBATCH --account=pet-gspt-2019
#SBATCH --partition=nodes
#SBATCH --time=02:00:00
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=200mb
ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Running GS2 code with ${SLURM_NTASKS} tasks and ${SLURM_CPUS_PER_TASK} cpus per task."
mpiexec -n ${SLURM_NTASKS} /users/mly509/gs2/bin/gs2 input.in>&out.txt
