import numpy as np

B=1.2


a = np.loadtxt(f"powerdata_{B}_0.txt", delimiter=",")

no_files = 1000

total_power = np.zeros([no_files,len(a)])

for ix in range(no_files):
  total_power[ix,:] = np.loadtxt(f"powerdata_{B}_{ix}.txt", delimiter=",")

power = sum(total_power)/no_files


import os
jobid = os.getenv('SLURM_ARRAY_TASK_ID')

np.savetxt(f'average_{B}.txt',power, delimiter=",")
