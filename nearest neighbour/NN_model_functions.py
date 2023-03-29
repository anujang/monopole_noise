# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 19:18:29 2021

@author: Anujan
"""


import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

# system at infinite temprature - max randomness
def state(N,M):
    spins = npr.randint(0,2, (N,M)) * 2 - 1  # produces number (-1,1) randomly
    return spins
    
# Dictionary of possible configurations
# s=1 is in ; s=-1 is out
# ground state
g_s1 = np.array([(1,1),
                 (-1,-1)]) ; g_s2 = g_s1*-1

g_s3 = np.array([(1,-1),
                 (1,-1)]) ; g_s4 = g_s3*-1

g_s5 = np.array([(1,-1),
                 (-1,1)]) ; g_s6 = g_s5*-1

# first excited state
e_s1 = np.array([(1,-1),
                 (1,1)]) ; e_s2 = e_s1*-1

e_s3 = np.array([(-1,1),
                 (1,1)]) ; e_s4 = e_s3*-1

e_s5 = np.array([(1,1),
                 (-1,1)]) ; e_s6 = e_s5*-1

e_s7 = np.array([(1,1),
                 (1,-1)]) ; e_s8 = e_s7*-1

e_s9 = np.array([(1,1),
                 (1,1)]) ; e_s10 = e_s9*-1

# to split larger matrices into 2x2 ones to identify each one from dictionary
def split(spins):
    n,m = np.shape(spins)
    num = int((n-(n%2))/2)  # works out how many 2x2 lattices in big lattice
    arr = np.arange(0,n,2)
    iso = np.zeros((num**2, 2, 2)) # empty array to put in 2x2 matrices, each around a lattice point
    ind=0
    for i in range(len(arr)):
        for j in range(len(arr)):
            indx = int(arr[i])
            indy = int(arr[j])
            iso[ind] = np.array([(spins[indx,indy], spins[indx,(indy+1)%n]), 
                          (spins[(indx+1)%n,indy], spins[(indx+1)%n,(indy+1)%n])])
            ind+=1
    return iso

def define(spins):
    gs1,gs2,gs3,gs4,gs5,gs6,es1,es2,es3,es4,es5,es6,es7,es8,es9,es10=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    for i in range(len(spins)):
        if np.all(spins[i] == g_s1):
            gs1+=1
        if np.all(spins[i] == g_s2):
            gs2+=1
        if np.all(spins[i] == g_s3):
            gs3+=1
        if np.all(spins[i] == g_s4):
            gs4+=1
        if np.all(spins[i] == g_s5):
            gs5+=1
        if np.all(spins[i] == g_s6):
            gs6+=1            
        if np.all(spins[i] == e_s1):
            es1+=1
        if np.all(spins[i] == e_s2):
            es2+=1
        if np.all(spins[i] == e_s3):
            es3+=1
        if np.all(spins[i] == e_s4):
            es4+=1
        if np.all(spins[i] == e_s5):
            es5+=1
        if np.all(spins[i] == e_s6):
            es6+=1
        if np.all(spins[i] == e_s7):
            es7+=1
        if np.all(spins[i] == e_s8):
            es8+=1
        if np.all(spins[i] == e_s9):
            es9+=1
        if np.all(spins[i] == e_s10):
            es10+=1
    return gs1,gs2,gs3,gs4,gs5,gs6,es1,es2,es3,es4,es5,es6,es7,es8,es9,es10

    
def energy(spins):
    c = spins
    E = 0
    #for i in range(len(c)):
    if c[0]>0 or c[1]>0 or c[2]>0 or c[3]>0:   # states gs1 to gs4 have energy -2J2
        E+=(c[0]+c[1]+c[2]+c[3])*((np.sqrt(2)-1)/(np.sqrt(2)-0.5))
        #print(E)
    if c[4]>0 or c[5]>0:                       # states gs5 and g6 have energy -4J1 + 2J2
        E+=(c[4]+c[5])*0#((-4*J1)+(2*J2))
        #print(E)
    if c[14]>0 or c[15]>0:                     # states es9 and es10 have energy 4J1 + 2J2
        E+=(c[14]+c[15])*((4*np.sqrt(2))/(2*np.sqrt(2) - 1))
        #print(E)
    if c[6]>0 or c[7]>0 or c[8]>0 or c[9]>0 or c[10]>0 or c[11]>0 or c[12]>0 or c[13]>0:
        E+=(c[6]+c[7]+c[8]+c[9]+c[10]+c[11]+c[12]+c[13]+c[14])*1
    return E
    
# import time
# timenow = time.time()
# a = state(2,2)
# b = split(a)
# c = define(b, len(a), len(a))
# d = energy(define,len(a), len(a))

# print(d)

#print("the energy in the system is", d, "J")
#timeend = time.time()
#print(timeend-timenow)