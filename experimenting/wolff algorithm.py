from copy import copy
import time 


timenow = time.time()

def getNN(site_indices, site_ranges, num_NN):
    '''
        site_indices: [i,j], site to get NN of
        site_ranges: [Nx,Ny], boundaries of the grid
        num_NN: number of nearest neighbors, usually 1

        function which gets NN on any d dimensional cubic grid
        with a periodic boundary condition
    '''

    Nearest_Neighbors = list();
    for i in range(len(site_indices)):
        for j in range(-num_NN,num_NN+1): #of nearest neighbors to include
            if(j == 0): continue;
            NN = list(copy(site_indices)); #don't want to overwite;
            NN[i] +=j;
            if(NN[i] >= site_ranges[i]):
                NN[i] = NN[i] - site_ranges[i];
            if(NN[i] < 0):
                NN[i] = site_ranges[i]+NN[i];
            Nearest_Neighbors.append(tuple(NN))
    return Nearest_Neighbors;


import numpy as np
from copy import copy
from scipy.ndimage.interpolation import shift

'''
    series of functions which will extract all major order parameters from the 
    Ising simulation
'''

def magnetization(lattice):
    N_s = np.prod(lattice.shape)
    return abs((1/N_s) * np.sum(lattice));

def magnetization_2(grid):
    ''' calculate <m^2>'''
    N_s = np.prod(grid.shape);
    return ((1/N_s))*np.sum(grid*grid); #if we do grid*grid, this is some constant number...

def energy(grid, J=1):
    '''
    this function must be able to generallize to arbitrary dimensions
    :param grid:
    :param J:
    :return:
    '''
    N = grid.shape;
    dimension = len(N); #dimension of system being simulated
    NN = copy(grid);
    E = 0;
    neighbors = 0;

    for i in range(dimension):
        for j in [-1,1]:
            neighbors += np.roll(NN, shift=j, axis = i);
            E+=J*np.sum(grid*np.roll(NN, shift=j, axis = i));
    DeltaE = J * (grid* neighbors)/(np.prod(N));
    return  -np.sum(DeltaE)/2;#-(E/np.prod(N))/2; #return is avg energy per site

def energy_2(grid,J=1):
    ''' measurement of energy^2 per site, which should be constant'''
    N_s = np.prod(grid.shape);
    N = grid.shape;
    dimension = len(N); #dimension of system being simulated
    NN = copy(grid);
    E = 0;
    for i in range(dimension):
        for j in [-1,1]:
            E+=J*(grid*np.roll(NN, shift=j, axis = i));
    return np.sum(E*E)/N_s/2;

def s0sr(grid, r):
    '''
    :param grid:
    :param r:
    :return:
    '''
    '''there is no special point so just pick one origin for o'''
    N_s = grid.shape;
    s0 = grid[0,0];
    ssr = 0;
    for i in range(N_s[0]):
        for j in range(N_s[1]):
            ssr += s0*grid[i,j];
    return ssr/np.prod(N_s)

def susceptibility(grid, beta):
    '''
    formula chi = 1/T(<m^2> - <m>^2)
    :param grid:
    :param beta:
    :return:
    '''
    m = magnetization(grid);
    m_2 = magnetization_2(grid);
    chi = beta*(m_2 - m**2);
    return chi;

def heat_capacity(grid,K):
    E = energy(grid, K);
    E2 = energy_2(grid, K);
    cv = K**2*(E2 - E**2);
    return cv;

import numpy as np
#from core_functions.getNN import getNN;
import matplotlib.pyplot as plt
#from variable_calculations.order_measures import *

## Wolf function
def Wolff_simulation(Lattice,K, epochs, thermalization_epochs = 100, num_views = 10):
    '''
    :param Lattice:
    :param K:
    :param epochs:
    :param thermalization_epochs:
    :param num_views:
    :return:
    '''
    plt.ion();
    ## wolff test using Frontier idea
    N = Lattice.shape;
    # generate random particle
    Lattice_History = list();
    p = 1 - np.exp(-2 * K)
    data = list();
    for t in range(epochs):
        change_tracker = np.ones(N);
        visited = np.zeros(N);
        root = []; # generate random coordinate by sampling from uniform random...
        for i in range(len(N)):
            root.append(np.random.randint(0, N[i], 1)[0])
        root = tuple(root);
        visited[root]=1;
        C = [root];  # denotes cluster coordinates
        F_old = [root];  # old frontier
        change_tracker[root] = -1;
        while (len(F_old) != 0):
            F_new = [];
            for site in F_old:
                site_spin = Lattice[tuple(site)]
                # get neighbors
                NN_list = getNN(site, N, num_NN=1);
                for NN_site in NN_list: ## if we do the full search, this is bad, because
                    nn = tuple(NN_site)
                    if (Lattice[nn] == site_spin and visited[nn] == 0):
                        if (np.random.rand() < p):
                            F_new.append(nn); visited[nn] = 1;
                            C.append(nn);
                            change_tracker[nn] = -1;
            F_old = F_new;

        # update the cluster
        Lattice = Lattice*change_tracker;
        if(t > thermalization_epochs):
            data.append(magnetization(Lattice));
        # for site in C:
        #     Lattice[site] = -1 * Lattice[site]
        # if (t % int(epochs/num_views) == 0):
        #     print('epoch: ' + str(t));
        #     plt.imshow(Lattice);
        #     plt.pause(0.05)

    # plt.imshow(Lattice);
    # plt.show()
    return Lattice, data;

def run_Wolff_epoch(Lattice, N, p):
    '''
    run one wolff epoch
    :param Lattice:
    :param N:
    :param p: 1-exp(-2K);
    :return:
    '''

    change_tracker = np.ones(N);
    visited = np.zeros(N);
    root = [];  # generate random coordinate by sampling from uniform random...
    for i in range(len(N)):
        root.append(np.random.randint(0, N[i], 1)[0])
    root = tuple(root);
    visited[root] = 1;
    C = [root];  # denotes cluster coordinates
    F_old = [root];  # old frontier
    change_tracker[root] = -1;
    while (len(F_old) != 0):
        F_new = [];
        for site in F_old:
            site_spin = Lattice[tuple(site)]
            # get neighbors
            NN_list = getNN(site, N, num_NN=1);
            for NN_site in NN_list:  ## if we do the full search, this is bad, because
                nn = tuple(NN_site)
                if (Lattice[nn] == site_spin and visited[nn] == 0):
                    if (np.random.rand() < p):
                        F_new.append(nn);
                        visited[nn] = 1;
                        C.append(nn);
                        change_tracker[nn] = -1;
        F_old = F_new;
    Lattice = Lattice * change_tracker;

    return Lattice;
## wolff test using Frontier idea
N = (512,512);
#lattice
Lattice = 2*np.random.randint(0,2,N)-1;
#Lattice = 2*np.ones(N);
#generate random particle
change_tracker = np.ones(N);

#K = 0.8;
p = 2/(np.log(1+np.sqrt(2)))#1-np.exp(-2*K)
epochs = 800; #epochs scales with number of sites that must be probed...
for t in range(epochs):
    root = [];
    for i in range(len(N)):
        root.append(np.random.randint(0,N[i],1)[0])
    root = tuple(root);
    #print('root: '+str(root))
    C = [root]; # denotes cluster coordinates
    visited = np.zeros(N);
    visited[root] = 1;
    F_old = [root]; #old frontier
    change_tracker[root] = -1;
    computations_tracker = 0;
    while(len(F_old) != 0): #there are situations where the cluster search slows down when the clusters get large
        F_new = [];
        for site in F_old:
            site_spin = Lattice[tuple(site)]
            #get neighbors
            NN_list = getNN(site, N, num_NN= 1);
            for NN_site in NN_list: ## this has to probe through 2*d nearest neighbors...
                computations_tracker+=1;
                nn = tuple(NN_site)
                if(Lattice[nn] == site_spin and visited[nn] == 0): # original codenn not in C):
                    ## nn not in C is an expensive computation if C is just a simple array
                    if(np.random.rand() < p):
                        F_new.append(nn);
                        C.append(nn); visited[nn] = 1;
                        change_tracker[nn] = -1;
        F_old = F_new;

    #at the end we have a cluster
    #print(C)
    #print(len(C))
    print('computations: '+str(computations_tracker))
    #update the cluster
    Lattice = Lattice*change_tracker;
    # for site in C:
    #     Lattice[site] = -1*Lattice[site]
    if(t%80 == 0):
        plt.imshow(Lattice);
        plt.show();

plt.imshow(Lattice, cmap='Greys');

plt.show()

timeend= time.time()

print(timeend - timenow)