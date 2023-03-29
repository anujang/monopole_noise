# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:56:55 2021

@author: Anujan
"""

import numpy as np
import numpy.random as npr
from numpy.random import rand
import matplotlib.pyplot as plt
import time
import scipy.constants as scc
import os
path = "C:\\Users\\Dell\\OneDrive - Cardiff University\\Year 4\\Project\\NN model"
os.chdir(path)
from NN_model_functions import state, split, define, energy

def initialstate(N):   
    ''' 
    Generates a random spin configuration for initial condition
    '''
    state = 2*np.random.randint(2, size=(N,N))-1
    return state

def mcmove(config, beta):
    '''
    Monte Carlo move using Metropolis algorithm 
    '''
    N=len(config)
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                s =  config[a, b]
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                cost = 2*s*nb
                
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                config[a, b] = s
    return config


init_state = state(2,2)
vertices = split(init_state)
configs = define(vertices)
ene = energy(configs)
print(init_state)
print(vertices)