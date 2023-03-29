import numpy as np


a = np.loadtxt(f"mags_0.txt", delimiter=",") 

no_files = 1000

total_power = np.zeros([no_files,len(a)]) 

for ix in range(no_files):
  total_power[ix,:] = np.loadtxt(f"mags_{ix}.txt", delimiter=",") 

mags = sum(total_power)/no_files 


import os 
jobid = os.getenv('SLURM_ARRAY_TASK_ID') 
np.savetxt(f'average_statics.txt', mags, delimiter=",")
