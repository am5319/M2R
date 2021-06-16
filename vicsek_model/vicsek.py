import numpy as np
from numpy.lib.function_base import average
import matplotlib.pyplot as plt

"""
A code to implement the vicsek model.
"""

# Set up all the parameters for the simulations.
L = 32
N = 500
R = 1
v = 0.03
eta = 0.1
delta_t = 1
final_t = 20
alpha = 0
rho = 4

# Pre-determined and random parameters
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))
R_sq = R**2
eta_list = np.arange(0, 5.1, 0.25)


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
                    if r_ijsq < R_sq and i != j:
                        neighbour_angles.append(directions[j])
        average_angle = average(neighbour_angles)
        average_angle_list.append(average_angle)
    average_angle_array = np.array(average_angle_list)
    directions = average_angle_array + np.random.uniform(-eta/2, eta/2, size=N)
    positions = np.mod(new_positions, L)
    return positions, directions


def simulate(positions, directions, L, N, R, delta_t, eta, v_list, final_t):
    for i in range(0, final_t, delta_t):
        positions, directions = step(positions, directions, L, N, R, delta_t,
                                     eta, v_list)
        v_x = v * np.cos(directions)
        v_y = v * np.sin(directions)
        v_list = np.column_stack((v_x, v_y))
    return positions, directions, v_x, v_y


def visualise(positions, directions, L, N, R, delta_t, eta, v_list, final_t):
    positions, directions, v_x, v_y = simulate(positions, directions, L, N, R,
                                               delta_t, eta, v_list, final_t)
    x_pos, y_pos = np.hsplit(positions, 2)

    plt.axis('equal')
    plt.quiver(x_pos, y_pos, v_x, v_y)
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    plt.show()


def order_param(positions, directions, L, N, R, delta_t, eta, v_list, final_t,
                v):
    new_positions, new_directions, new_v_x, new_v_y = simulate(
        positions, directions, L, N, R, delta_t, eta, v_list, final_t)
    sum_new_v_x = sum(new_v_x)
    sum_new_v_y = sum(new_v_y)
    order_param = (1 / (N * v)) * np.sqrt((sum_new_v_x ** 2) +
                                          (sum_new_v_y ** 2))
    return order_param


def step_v2(positions, directions, L, N, R, delta_t, eta, v_list, alpha):
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
                    if r_ijsq < R_sq and i != j:
                        pos_ij = positions[j] - positions[i]
                        cos_ang = np.dot(pos_ij, v_list[i])
                        weight = (1 + (alpha * np.sign(cos_ang))) / 2
                        neighbour_angles.append(weight * directions[j])
        average_angle = average(neighbour_angles)
        average_angle_list.append(average_angle)
    average_angle_array = np.array(average_angle_list)
    directions = average_angle_array + np.random.uniform(-eta/2, eta/2, size=N)
    positions = np.mod(new_positions, L)
    return positions, directions


def simulate_v2(positions, directions, L, N, R, delta_t, eta, v_list, final_t,
                alpha):
    for i in range(0, final_t, delta_t):
        positions, directions = step_v2(positions, directions, L, N, R,
                                        delta_t, eta, v_list, alpha)
        v_x = v * np.cos(directions)
        v_y = v * np.sin(directions)
        v_list = np.column_stack((v_x, v_y))
    return positions, directions, v_x, v_y


def visualise_v2(positions, directions, L, N, R, delta_t, eta, v_list, final_t,
                 alpha):
    positions, directions, v_x, v_y = simulate_v2(positions, directions, L, N,
                                                  R, delta_t, eta, v_list,
                                                  final_t, alpha)
    x_pos, y_pos = np.hsplit(positions, 2)

    plt.axis('equal')
    plt.quiver(x_pos, y_pos, v_x, v_y)
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    plt.show()
