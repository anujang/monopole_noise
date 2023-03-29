import numpy as np
from numpy.random import rand
import os
path = "C:\\Users\\Dell\\OneDrive - Cardiff University\\Year 4\\Project\\NN model"
os.chdir(path)
from NN_model_functions import state, split, define, energy

def calcMag(config):
    '''
    Magnetization of a given configuration
    '''
    mag = np.sum(config)
    return mag


init_state = state(14,14)
print("before",calcMag(init_state))
beta = 100
N=len(init_state)
for i in range(N):
    for j in range(N):
            a = np.random.randint(0, N)
            b = np.random.randint(0, N)
            s =  init_state[a, b]
            
            if a%2 == 0 and b%2 == 0:     # top left of 2x2
                TL = np.array([[init_state[a,b],init_state[a,b+1]],[init_state[a+1,b],init_state[a+1,b+1]]]) # the 2x2 array within which the chosen point resides
                
                # if the TL point is chosen, the corresponding vertex where the chosen point is its BR point
                BR_c = np.array([[init_state[(a-1)%2,(b-1)%2],init_state[(a-1)%2,b]],[init_state[a,(b-1)%2],init_state[a,b]]])
                energy_before = energy(define(split(TL))) + energy(define(split(BR_c)))
                print(energy_before)
                nb = init_state[(a+1)%N,b] + init_state[a,(b+1)%N] + init_state[(a-1)%N,b] + init_state[a,(b-1)%N]
                cost = 2*s*nb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                TL[0,0] = s
                BR_c[1,1] = s
                
                energy_after = energy(define(split(TL))) + energy(define(split(BR_c)))
                #print("change in energy:", energy_after - energy_before)

                #print("before",calcMag(init_state))
                init_state[a,b] = TL[0,0]
                init_state[a,b+1] = TL[0,1]
                init_state[a+1,b] = TL[1,0]
                init_state[a+1,b+1] = TL[1,1]
                init_state[(a-1)%2,(b-1)%2] = BR_c[0,0]
                init_state[(a-1)%2,b] = BR_c[0,1]
                init_state[a,(b-1)%2] = BR_c[1,0]
                init_state[a,b] = BR_c[1,1]
                #print("after",calcMag(init_state))

                
            elif a%2 == 0 and b%2 == 1:   # top right of 2x2
                TR = np.array([[init_state[a,b-1],init_state[a,b]],[init_state[a+1,b-1],init_state[a+1,b]]])
                
                # if the TR point is chosen, the corresponding vertex where the chosen point is its BL point
                BL_c = np.array([[init_state[(a-1)%2,b],init_state[(a-1)%2,(b+1)%2]],[init_state[a,b],init_state[a,(b+1)%2]]])
                energy_before = energy(define(split(TR))) + energy(define(split(BL_c)))
                
                nb = init_state[(a+1)%N,b] + init_state[a,(b+1)%N] + init_state[(a-1)%N,b] + init_state[a,(b-1)%N]
                cost = 2*s*nb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                TR[0,1] = s
                BL_c[1,0] = s
                
                energy_after = energy(define(split(TR))) + energy(define(split(BL_c)))
                #print("change in energy:", energy_after - energy_before)

                #print("before",calcMag(init_state))
                init_state[a,b-1] = TR[0,0]
                init_state[a,b] = TR[0,1]
                init_state[a+1,b-1] = TR[1,0]
                init_state[a+1,b] = TR[1,1]
                
                init_state[(a-1)%2,b] = BL_c[0,0]
                init_state[(a-1)%2,(b+1)%2] = BL_c[0,1]
                init_state[a,b] = BL_c[1,0]
                init_state[a,(b+1)%2] = BL_c[1,1]
                #print("after",calcMag(init_state))

                
                


            elif a%2 == 1 and b%2 == 0:   # bottom left of 2x2
                BL = np.array([[init_state[a-1,b],init_state[a-1,b+1]],[init_state[a,b],init_state[a,b+1]]])
                
                # if the BL point is chosen, the corresponding vertex where the chosen point is its TR point
                TR_c = np.array([[init_state[a,(b-1)%2],init_state[a,b]],[init_state[(a+1)%2,(b-1)%2],init_state[(a+1)%2,b]]])
                energy_before = energy(define(split(BL))) + energy(define(split(TR_c)))
                
                nb = init_state[(a+1)%N,b] + init_state[a,(b+1)%N] + init_state[(a-1)%N,b] + init_state[a,(b-1)%N]
                cost = 2*s*nb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                BL[0,1] = s
                TR_c[1,0] = s
                
                energy_after = energy(define(split(BL))) + energy(define(split(TR_c)))
                #print("change in energy:", energy_after - energy_before)
                
               
                #print("before",calcMag(init_state))
                init_state[a-1,b] = BL[0,0]
                init_state[a-1,b+1] = BL[0,1]
                init_state[a,b] = BL[1,0]
                init_state[a,b+1] = BL[1,1]
                
                init_state[a,(b-1)%2] = TR_c[0,0]
                init_state[a,b]  = TR_c[0,1]
                init_state[(a+1)%2,(b-1)%2] = TR_c[1,0]
                init_state[(a+1)%2,b] = TR_c[1,1]
               # print("after",calcMag(init_state))
               
                
            elif a%2 == 1 and b%2 == 1:   # bottom right of 2x2
                BR = np.array([[init_state[a-1,b-1],init_state[a-1,b]],[init_state[a,b-1],init_state[a,b]]])
                
                # if the BR point is chosen, the corresponding vertex where the chosen point is its TL point
                TL_c =  np.array([[init_state[a,b],init_state[a,(b+1)%2]],[init_state[(a+1)%2,b],init_state[(a+1)%2,(b+1)%2]]])
                energy_before = energy(define(split(BR))) + energy(define(split(TL_c)))
                
                nb = init_state[(a+1)%N,b] + init_state[a,(b+1)%N] + init_state[(a-1)%N,b] + init_state[a,(b-1)%N]
                cost = 2*s*nb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                BR[0,1] = s
                TL_c[1,0] = s
                
                energy_after = energy(define(split(BR))) + energy(define(split(TL_c)))
                #print("change in energy:", energy_after - energy_before)
                
                
               # print("before",calcMag(init_state))
                init_state[a-1,b-1] = BR[0,0]
                init_state[a-1,b] = BR[0,1]
                init_state[a,b-1] = BR[1,0]
                init_state[a,b] = BR[1,1]
                
                init_state[a,b] = TL_c[0,0]
                init_state[a,(b+1)%2] = TL_c[0,1]
                init_state[(a+1)%2,b] = TL_c[1,0]
                init_state[(a+1)%2,(b+1)%2] = TL_c[1,1]
               # print("after",calcMag(init_state))
print("after",calcMag(init_state))