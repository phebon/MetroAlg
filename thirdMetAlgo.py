import numpy as np
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import datetime as dt
import pickle
import math
import data_handling
#import thread
start_time = time.time()

""" This initializes a initial cell by either assigning each point on the grid
     a random element of {-1,1} (uniform distribution) or a homogenous Matrix of 1s. 
     This behaviour is controlled by the <rand> argument. 
     The resulting matrix is padded by two columns and rows of 0s to each side, to simplify the algorithm. 
     We therefore assume, that we have non periodic boundary conditions."""
def initial_cell(n,rand=False):
    if rand:
        bin_matr = np.random.randint(0,2,(n,n))
        sign_matr = list(map(lambda x: 2*x-1,bin_matr))
        return np.pad(sign_matr,(2,2),"constant",constant_values = (0,0))
    else:
        return np.pad(np.ones((n,n)),(1,1),"constant",constant_values = (0,0))
# This function gets all the neighbours of a given position <pos>
def neighbour_indices(pos):
    x,y = pos
    return ((x+1,y),(x-1,y),(x,y-1),(x,y+1))

"""This calculates the Hamiltonian of one gridpoint by adding up the interactions
    with the 4 closest neighbours of the cell, 
    which simplify to 2 or three neighbours for the edge and corner cases respectively."""
def single_ham(pos,matrix):
    neighbour_list = neighbour_indices(pos)
    #print('neighbours are: ', neighbour_list)
    #here we can speed up by indexing with a single index.
    #products = [matrix[el[0]][el[1]]*matrix[pos[0]][pos[1]] for el in neighbour_list]
    #print('products are: ',products)
    ham_pos = np.sum([matrix[el[0]][el[1]]*matrix[pos[0]][pos[1]] for el in neighbour_list])
    #print("ham_pos is ",ham_pos,f"for {pos}")
    return ham_pos

"""This is the relevant part of the hamiltonian of a grid, 
    which matters for the difference of Energy between two states. 
    It is two times the hamiltoinian of the cell, that is supposed to be changed. """
"""def change_ham(pos, matrix):
    nbrs = neighbour_indices(pos)
    diff_energy = np.sum([single_ham(el, matrix) for el in nbrs])
    return diff_energy"""

def change_ham(pos,matrix):
    nbrs = neighbour_indices(pos)
    diff_energy = 2*single_ham(pos, matrix)
    return diff_energy

def act_on(pos, matrix, beta):
    x,y = pos
    matrix_copy = np.copy(matrix)
    matrix_copy[x][y]=-matrix[x][y]
    energy_diff = 2*single_ham((x,y),matrix)#-change_ham((x,y),matrix)
    #print("Energy Difference is",energy_diff)
    #print("first: \n",matrix,"\nThen:\n", matrix_copy)
    if energy_diff < 0:
        #print('flipped!')
        '''print("orig. matrix:\n",matrix,"\nnew matrix:\n",
                                            matrix_copy,
                                            f"\nchanged element: {matrix[x][y]}->{matrix_copy[x][y]} at {x},{y}")
                                contin = input("Nächster Zug?")'''
        return matrix_copy
    else:
        xrand = np.random.uniform()
        if xrand<=math.e**(- beta * energy_diff):
            '''print("orig. matrix:\n",matrix,"\nnew matrix:\n",
                                                            matrix_copy,
                                                            f"\nchanged element: {matrix[x][y]}->{matrix_copy[x][y]} at {x},{y}, xrand = {xrand}<{np.exp(-beta*energy_diff)}")
                                                contin = input("Nächster Zug?")'''
            return matrix_copy
        else:
            #print(f"xrand = {xrand} > {np.exp(-beta*energy_diff)}")
            return matrix

    #"""elif np.random.uniform() <= np.exp(- beta * energy_diff):
     #               #print('flipped!')
      #              print("orig. matrix:\n",matrix,"\nnew matrix:\n",
       #                 matrix_copy,
        #                f"\nchanged element: {matrix[x][y]}->{matrix_copy[x][y]} at {x},{y}")
         #           contin = input("Nächster Zug?")
          #          return matrix_copy
           #     else:
            #        #print(np.exp(-energy_diff))
             #       return matrix"""
def matrix_update(matrix, n, beta):
    return act_on( np.random.randint(1,n+1,2), matrix, beta)


def sweeps(matrix, n, beta):
    fir = matrix
    for ind in range(n*n):
        sec = matrix_update(fir, n, beta)
        fir = sec
    return sec


# This is a function to thermalize a nxn matrix.
def thermalize_matrix(matrix, n, beta):
    ind = 0
    fir = matrix
    for ind in range(12*n**2):
        if ind%500 ==0: print(f"thermosweep number: {ind}")
        sec = sweeps(fir, n, beta)
        fir = sec
    return sec


def get_observable_value(func, matrix, n, beta):
    return func(matrix, n, beta)


def measure_obs(funcs, matrix, beta, n, number_of_measurements):
    thermo = thermalize_matrix(matrix, n, beta)
    obs_val_list = [[func(matrix, n, beta) for func in funcs]]
    curr_matrix = thermo
    for ind in range(number_of_measurements-1):
        if ind%500 ==0: print(f"Measurement number {ind}")
        for inde in range(10):
            new_matrix = sweeps(curr_matrix, n, beta)
            curr_matrix = new_matrix
        obs_val_list.append([func(curr_matrix, n, beta) for func in funcs])
    #print(curr_matrix)
    return obs_val_list


# The magnetization is normalized to n^2, which is the particle number.
def magnetization(matrix, n, beta):
    return sum(matrix.flatten())/(n**2)

def susceptibility_mom(matrix, n, beta):
    mag = magnetization(matrix, n, beta)
    sus_mom = 0
    for xindex in range(1,n+1):
        for yindex in range(1,n+1):
            sus_mom += beta*(matrix[xindex][yindex]-mag)**2
    return sus_mom 

def average(arr):
    return sum(arr)/len(arr)

'''Here we set up our measurements'''
n           = 20
num_of_meas = 4000
#betas    = list(reversed([0.1,0.3,0.5]+[.6+i*.05 for i in range(8)]+[.98,.99]))
#betas = [0.15+0.05*i for i in range(11)]
#betas = [0.4+0.01*i for i in range(16)]#betas    =[0,0.5,1,2,3,4]
#betas = [.1,.2,.3,.32,.34,.38,.4,.42,.44,.48,.5,.6,.7,.8,.9]
betas = [.1,.15,.2,.3,.32,.34,.36,.38,.4,.42,.43,.44,.46,.48,.5,.55,.6,.65,.7,.8,.9]
betas = [.25,.28,.41,.45,.47,.49,.52,.75,.85]
#betas    = list(reversed([0.1]))
#betas    = [0.1,0.2,0.4,0.5]
exp_list = []
#thread.start_new_thread( print_time, ("Thread-1", 2, ) )
#thread.start_new_thread( measure_obs(), ("Thread-2", 4, ) )
for beta_el in betas:
    loop_start = time.time()
    print("beta = ", beta_el)

    obs_list = measure_obs([magnetization, susceptibility_mom], initial_cell(n), beta_el, n, num_of_meas)    

    avg_magn    = [np.mean([abs(el[0]) for el in obs_list]),np.std([el[0] for el in obs_list])/(num_of_meas**(.5))]
    avg_sus     = [np.mean([el[1] for el in obs_list]),np.std([el[1] for el in obs_list])/(num_of_meas**.5)]

    exp_list.append([avg_magn, avg_sus])
    print(f"This took {loop_start-time.time()} seconds for thermalizing and {num_of_meas} measurements!")

fig, ax = plt.subplots(2,figsize=(10,10))
# Here come the magnetization and susceptibility lists:
magn_list = [el[0][0] for el in exp_list]
magn_sig_list = [el[0][1] for el in exp_list]

sus_list = [el[1][0] for el in exp_list]
sus_sig_list = [el[1][1] for el in exp_list]

file_name = data_handling.make_file_name(n, num_of_meas)
data_handling.save_data(file_name, betas, magn_list, magn_sig_list, sus_list, sus_sig_list)

print(magn_list)
ax[0].errorbar(betas, magn_list,magn_sig_list, fmt='x', color='b')
ax[1].errorbar(betas, sus_list,sus_sig_list, fmt='x', color='r')
fig.savefig(f'Results/plot_' + file_name + '.png', dpi=300)
print("<M>= ", avg_magn)
end_time = time.time()
print(f"The runtime of this script is {end_time-start_time} seconds")