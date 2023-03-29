#!/bin/bash --login
###
# job name
#SBATCH -J monopole_test

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --mail-user=anujan0902@gmail.com


#SBATCH --qos=maxjobs1500

#SBATCH --partition=htc

#SBATCH --array=0-999


#SBATCH -n 1

# maximum job time in D-HH:MM
#SBATCH --time=0-00:20

# number of tasks you are requesting
#SBATCH --ntasks=1

# memory per process in MB
#SBATCH --mem=1000

# number of nodes needed
#SBATCH --nodes=1

# specify our current project
# change this for your own work
#SBATCH --account=scw1887

#echo SLURM_JOB_ID $SLURM_JOB_ID
#echo SLURM_ARRAY_JOB_ID $SLURM_ARRAY_JOB_ID
#echo SLRUM_ARRAY_TASK_ID $SLRUM_ARRAY_TASK_ID

module load anaconda/2019.03
source activate scw_test
python3 monopole_single.py


