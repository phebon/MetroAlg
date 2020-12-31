import numpy as np
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import datetime as dt
import pickle
import math
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
        return np.pad(np.ones((n,n)),(2,2),"constant",constant_values = (0,0))
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
    energy_diff = change_ham((x,y),matrix_copy)-change_ham((x,y),matrix)
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
    return act_on( np.random.randint(2,n-1,2), matrix, beta)


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
        sec = sweeps(fir, n, beta)
        fir = sec
    return sec

def get_observable_value(func, matrix, n):
    return func(matrix, n)


def measure_obs(funcs, matrix, beta, n, number_of_measurements):
    thermo = thermalize_matrix(matrix, n, beta)
    obs_val_list = [[func(matrix, n) for func in funcs]]
    curr_matrix = thermo
    for ind in range(number_of_measurements-1):
        new_matrix = sweeps(curr_matrix, n, beta)
        curr_matrix = new_matrix
        obs_val_list.append([func(curr_matrix, n) for func in funcs])
    print(curr_matrix)
    return obs_val_list

def magnetization(matrix, n):
    return abs(sum(matrix.flatten()))

def average(arr):
    return sum(arr)/len(arr)

'''Here we set up our measurements'''
n       = 10
num_of_meas = 2000
betas    = list(reversed([0.1,0.3,0.5]+[.6+i*.05 for i in range(8)]+[.98,.99]))
betas = [0.15+0.025*i for i in range(11)]
#betas    =[0,0.5,1,2,3,4]
#betas    = list(reversed([0.1]))
mags_list = []
for beta_el in betas:
    loop_start = time.time()
    print("beta = ",{beta_el} )
    mag_list = measure_obs([magnetization], initial_cell(n), beta_el, n, num_of_meas)
    print("length of the mag_lsit is:", len(mag_list)," \nthe mag_list is:\n", mag_list[1000:1010])
    avg_magn = [np.mean([el[0] for el in mag_list]),np.std([el[0] for el in mag_list])]
    mags_list.append(avg_magn)
    print(f"This took {loop_start-time.time()} seconds for thermalizing and {num_of_meas} measurements!")
fig, ax = plt.subplots(figsize=(10,10))
magn_list = [el[0] for el in mags_list]
magn_sig_list = [el[1] for el in mags_list]
print(magn_list)
ax.errorbar(betas, magn_list,magn_sig_list, fmt='x', color='b')
#ax.set_xscale("log")
fig.savefig('Plots/first_try03.png', dpi=300)
#fir = initial_cell(10)
#sec = thermalize_matrix(fir,n,.5)
#third = sweeps(sec,n,.1)
#print(fir,sec,third,sec==third)
print("<M>= ", avg_magn)
end_time = time.time()
"""fir = initial_cell(6,True)
for ind in range(10):
    sec = act_on((3,3),fir,0.5,)
    fir = sec
print(f"The runtime of this script is {end_time-start_time} seconds")"""