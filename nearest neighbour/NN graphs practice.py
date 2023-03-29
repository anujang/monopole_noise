# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 22:32:16 2021

@author: Anujan
"""

import os
import numpy as np
import matplotlib.pyplot as plt
path = "C:\\Users\\Dell\\OneDrive - Cardiff University\\Year 4\\Project\\NN model"
os.chdir(path)
from NN_model_functions import state, split, define, energy
from rajeshrinet import initialstate, mcmove, calcEnergy, calcMag
import time
import scipy.constants as scc

nt = 32
L = 10
eqSteps= 2**4
mcSteps = 2**5


a_w = np.linspace(1.2, 2.4, nt);
e1 = ((np.sqrt(2)-1)/(np.sqrt(2)-0.5))
T = (-e1/(np.log(a_w)))
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)

g = np.zeros((nt,mcSteps))

ti = time.time()

for tt in range(nt):
    a = state(L,L)
    E1 = M1 = E2 = M2 = 0
    iT=1.0/T[tt]; iT2=iT*iT;
    print(tt)
    for i in range(eqSteps):
        mcmove(a, iT)
        
    for i in range(mcSteps):
        mcmove(a, iT)
        b = split(a)
        c = define(b, len(a), len(a))
        d = energy(define,len(a), len(a))
        g[tt,i] = np.abs(np.sum(a)/(L**2))
        #print(d)
        #print("the energy in the system is", d, "J")
        

# ti_end = time.time()
# print(ti_end - ti, "seconds")

# plt.figure()
# plt.plot(a_w, g[:,-1])

















# x = np.zeros(100)
# for ix in range(len(x)):
        
#     a = state(6,6)
#     b = split(a)
#     c = define(b, len(a), len(a))
#     d = energy(define,len(a), len(a))
#     x[ix] = calcMag(a)/(len(a)**2)
#     #print("the energy in the system is", d, "J")

# plt.plot(x)



