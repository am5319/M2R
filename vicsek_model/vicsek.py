import numpy as np
from numpy.lib.function_base import average
import matplotlib.pyplot as plt

"""
A code to implement the vicsek model.
"""

# Set up all the parameters for the simulations.
L = 10
N = 100
R = 1
v = 1
delta_t = 1
eta = 0.2
final_t = 20

# Pre-determined and random parameters
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))
R_sq = R**2


def step(positions, directions, L, N, R, delta_t, eta, v_list):
    average_angle_list = []
    new_positions = positions + (delta_t * v_list)
    for i in range(N):
        neighbour_angles = []
        for j in range(N):
            x_ij = positions[i][0] - positions[j][0]
            if abs(x_ij) < R:
                y_ij = positions[i][1] - positions[j][1]
                if abs(y_ij) < R:
                    r_ijsq = x_ij**2 + y_ij**2
                    if r_ijsq < R_sq:
                        neighbour_angles.append(directions[j])
        average_angle = average(neighbour_angles)
        average_angle_list.append(average_angle)
    average_angle_array = np.array(average_angle_list)
    directions = average_angle_array + np.random.uniform(-eta/2, eta/2, size=N)
    positions = new_positions % L
    return positions, directions


def simulate(positions, directions, L, N, R, delta_t, eta, v_list, final_t):
    for i in range(0, final_t, delta_t):
        positions, directions = step(positions, directions, L, N, R, delta_t,
                                     eta, v_list)
        v_x = v * np.cos(directions)
        v_y = v * np.sin(directions)
        v_list = np.column_stack((v_x, v_y))
    return positions, directions, v_x, v_y


positions, directions, v_x, v_y = simulate(positions, directions, L, N, R,
                                           delta_t, eta, v_list, final_t)
x_pos, y_pos = np.hsplit(positions, 2)

plt.quiver(x_pos, y_pos, v_x, v_y)
