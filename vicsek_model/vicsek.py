import numpy as np
from numpy.lib.function_base import average

"""
A code to implement the vicsek model.
"""

#Set up all the parameters for the simulations
L = 10
N = 100
r = 1
v = 1
delta_t = 1
eta = 0.2
final_t = 20

positions = np.random.uniform(0, L, size=(N, 2))
theta_list = np.random.uniform(-np.pi, np.pi, size=N)
t0 = 0
v_x = v * np.cos(theta_list)
v_y = v * np.sin(theta_list)
v_list = np.column_stack((v_x, v_y))

def step(N, r, delta_t, eta, positions, theta_list):
    alignment_terms = []
    new_positions = positions + (delta_t * v_list)
    for i in range(N):
        neighbour_angles = []
        for j in range(N):
            if ((positions[i][0] - positions[j][0])**2 + (positions[i][1] - positions[j][1])**2 <= r**2) and (i != j):
                neighbour_angles.append(theta_list[j])
        average_neighbour = np.average(neighbour_angles)
        theta_list[i] = average_neighbour + np.random.uniform(-eta/2, eta/2)
    positions = new_positions
    return positions
    return theta_list
