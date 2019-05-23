# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 12:30:17 2018

@author: Kyle
"""
import numpy as np
from read_adf11 import read_adf11
import matplotlib.pyplot as plt
from math import exp

n_states = 3
atom = "He"
dens, temps, acd_dat = read_adf11('acd96_he.dat')
dens, temps, scd_dat = read_adf11('scd96_he.dat')

n_temps = len(temps)
states = []
for i in range(n_states):
    states.append(atom+"{}".format(i))
    if (i != 0):
        states[i] += "+"

n_reps = 1000  # number of random samples of rates to take


def generate_matrix(temp_index, dens_index, scd_dat, acd_dat):
    matr = np.zeros((n_states, n_states))
    R = acd_dat[:,temp_index, dens_index]
    S = scd_dat[:, temp_index, dens_index]
    index = np.arange(n_states-1)
    
    matr[index, index] = -S[index]-R[index-1]
    matr[index+1, index] = S[index]
    matr[index, index+1] = R[index]
    matr[0,0] = -S[0]
    matr[n_states-1, n_states-1] = -R[n_states-2]
    matr *= dens[dens_index]
    return matr
    
    
def rand_matrix(matrix, frac):
    delta = matrix
    delta *= frac
    rand = 2*np.random.random(matrix.shape) - 1
    delta *= rand
    matrix += delta
    return matrix
    

def steady_state(matrix):
    matrix[0,:] = 1
    inverse = np.linalg.inv(matrix)
    vec = np.zeros(n_states)
    vec[0] = 1.0
    return np.dot(inverse, vec)

def time_evol_euler(matrix, t_final, t_steps):
    delta_t = t_final / t_steps
    time = np.linspace(0, t_final, t_steps)
    N = np.zeros((t_steps, n_states))
    N[0, 0] = 1
    for i in range(1, t_steps):
        if i==2:
            print(np.dot(matrix, N[i-1,:])*delta_t)
        N[i, :] = N[i-1, :] + np.dot(matrix, N[i-1, :]) * delta_t
    return (time, N)

def exact_soln(eig, coeff, initial, t):
    eig_vals = eig[0]
    eig_vecs = eig[1]
    solution = np.zeros(n_states)
    for i in range(n_states):
        solution += coeff[i] * eig_vecs[:,i] * exp(eig_vals[i]*t)
    return solution

state_vectors = np.zeros((n_temps, n_states))
for i in range(n_temps):
    state_vectors[i,:] = steady_state(generate_matrix(i,0, scd_dat, acd_dat))
plt.plot(temps, state_vectors)
plt.xscale("log")
plt.xlim(1, 100)

t_final = 100
t_steps = 1000
time, N = time_evol_euler(generate_matrix(27, 0, scd_dat, acd_dat), t_final, t_steps)
plt.figure()
for i in range(n_states):
    plt.semilogx(time, N[:,i])

print(generate_matrix(27, 0, scd_dat, acd_dat))
"""
frac = 0.1
t_final = 10000000000
t_steps = 1000
time_arr = np.linspace(0, t_final, t_steps)

N_arr = np.zeros((n_reps, t_steps, n_states))


for i in range(n_reps):
   
    initial = np.zeros(n_states)
    initial[0] = 1
    scd = rand_matrix(scd_dat, frac)
    acd = rand_matrix(acd_dat, frac)
    eig = np.linalg.eig(matrix(29, 20, scd, acd))
    coeff = np.linalg.solve(eig[1], initial)

    for index, t in enumerate(time_arr):
        N_arr[i, index,:] = exact_soln(eig, coeff, initial, t)
        

avg = np.mean(N_arr, axis=0)
err = np.std(N_arr, axis=0) /avg

fig, axes = plt.subplots(2, 1)
axes[0].semilogx(time_arr, avg)
axes[1].semilogx(time_arr, err)
axes[0].legend(states)
axes[1].legend(states)
"""
"""
def exact_soln(matrix, t):
    eig = np.linalg.eig(matrix)
    eig_vals = eig[0]
    eig_vecs = eig[1]
    initial = np.zeros(n_states)
    initial[0] = 1
    coeff = np.linalg.solve(eig_vecs, initial)
    solution = np.zeros(n_states)
    for i in range(n_states):
        solution += coeff[i] * eig_vecs[:,i] * exp(eig_vals[i]*t)
    return solution
"""



"""
fig, axes = plt.subplots(2, 1)

for i in range(n_states):
    axes[0].semilogx(temps, state_vectors[:, i])
for i in range(n_states):
    axes[1].semilogx(time_arr, N_arr[:, i])
"""

