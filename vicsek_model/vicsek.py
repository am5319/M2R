import numpy as np

"""
A code to implement the vicsek model for flocks of birds, schools of fish, etc.
"""

#Set up all the parameters for the simulations
L = 10
N = 100
r0 = 1
delta_t = 1
final_t = 20
eta = 1

positions = np.random.uniform(0, L, size=(N, 2))
theta_list = np.random.uniform(0, 2*np.pi, size=N)
