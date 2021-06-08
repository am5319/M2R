import numpy as np
from numpy.lib.function_base import average
import matplotlib.pyplot as plt


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


# Set up all the parameters for the simulations.
rho = 4
R = 1
v = 0.03
delta_t = 1
final_t = 20

# Pre-determined and random parameters
R_sq = R**2


eta_list = np.arange(0, 5, 0.1)

N, L = 50, np.sqrt(50 / rho)
order_param_list_50 = []
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))

for i in eta_list:
    new_positions, new_directions, new_v_x, new_v_y = simulate(
        positions, directions, np.sqrt(12.5), 50, R, delta_t, i, v_list,
        final_t)
    new_v_x_sum = sum(new_v_x)
    new_v_y_sum = sum(new_v_y)
    order_param = (1 / (50 * v)) * np.sqrt(new_v_x_sum**2 + new_v_y_sum**2)
    order_param_list_50.append(order_param)


N, L = 100, np.sqrt(100 / rho)
order_param_list_100 = []
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))

for i in eta_list:
    new_positions, new_directions, new_v_x, new_v_y = simulate(
        positions, directions, 5, 100, R, delta_t, i, v_list, final_t)
    new_v_x_sum = sum(new_v_x)
    new_v_y_sum = sum(new_v_y)
    order_param = (1 / (100 * v)) * np.sqrt(new_v_x_sum**2 + new_v_y_sum**2)
    order_param_list_100.append(order_param)


N, L = 500, np.sqrt(500 / rho)
order_param_list_500 = []
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))

for i in eta_list:
    new_positions, new_directions, new_v_x, new_v_y = simulate(
        positions, directions, np.sqrt(125), 500, R, delta_t, i, v_list,
        final_t)
    new_v_x_sum = sum(new_v_x)
    new_v_y_sum = sum(new_v_y)
    order_param = (1 / (500 * v)) * np.sqrt(new_v_x_sum**2 + new_v_y_sum**2)
    order_param_list_500.append(order_param)


N, L = 5000, np.sqrt(5000 / rho)
order_param_list_5000 = []
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))

for i in eta_list:
    new_positions, new_directions, new_v_x, new_v_y = simulate(
        positions, directions, np.sqrt(1250), 5000, R, delta_t, i, v_list,
        final_t)
    new_v_x_sum = sum(new_v_x)
    new_v_y_sum = sum(new_v_y)
    order_param = (1 / (5000 * v)) * np.sqrt(new_v_x_sum**2 + new_v_y_sum**2)
    order_param_list_5000.append(order_param)


N, L = 10000, np.sqrt(10000 / rho)
order_param_list_10000 = []
positions = np.random.uniform(0, L, size=(N, 2))
directions = np.random.uniform(-np.pi, np.pi, size=N)
v_x = v * np.cos(directions)
v_y = v * np.sin(directions)
v_list = np.column_stack((v_x, v_y))

for i in eta_list:
    new_positions, new_directions, new_v_x, new_v_y = simulate(
        positions, directions, 50, 10000, R, delta_t, i, v_list, final_t)
    new_v_x_sum = sum(new_v_x)
    new_v_y_sum = sum(new_v_y)
    order_param = (1 / (10000 * v)) * np.sqrt(new_v_x_sum**2 + new_v_y_sum**2)
    order_param_list_10000.append(order_param)


plt.plot(eta_list, order_param_list_50, "s")
plt.plot(eta_list, order_param_list_100, "+")
plt.plot(eta_list, order_param_list_500, "x")
plt.plot(eta_list, order_param_list_5000, "^")
plt.plot(eta_list, order_param_list_10000, "D")
plt.show()
