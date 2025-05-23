"""
1D Rocket Flight Simulator
==========================

Author: Ivan Krat 
Initial Release: 2021  
Latest Update: 2025-05-22  
Version: 3.2

Description:
------------
This script simulates the vertical (1D) flight of a rocket, accounting for:
- Gravity
- Aerodynamic drag
- Engine burn time and thrust

The simulator tracks the rocket's altitude, velocity, and acceleration over time,
providing a simplified but realistic model of suborbital rocket flight in Earth's atmosphere.

Libraries Used:
---------------
- NumPy: for efficient numerical computation
- Matplotlib: for plotting flight trajectory

Assumptions:
------------
- Atmospheric density is treated as a function of altitude.
- Drag force is calculated using a fixed or parameterized drag coefficient.
- Thrust is assumed to be constant during the burn period.
- The rocket is constrained to vertical motion (no pitch/yaw dynamics).

Version History:
----------------
v1.0 - 2021: Initial release with basic motion under constant thrust  
v2.0 - Added aerodynamic drag modeling  
v3.0 - Replaced loop-based calculations with a more efficient vectorized approach using linear algebra
"""



import matplotlib.pyplot as plt
import numpy as np



# Constants:
rho_STP = 23.67*(10**-4)  # (desnity at standard tempurature and pressure) slug/ft^3
g = 32.17                 # (gravity) ft/s^2

# Input Variables
A = 10**2/4*np.pi/144     # (frontal area) ft^2
Cd = 0.75                 # (ballistic coeff)
F_t = 5000                # (thrust) lbm
m_0 = 200/g               # (dry mass) slugs
m_p = 250/g               # (propellent mass) slugs
m_dot = 15/g              # (mass flow) slugs
t_step = 0.1              # (delta t) sec

# Calculated Values
t_burn = m_p/m_dot        # (burn time) sec
m_net = m_0 + m_p         # (starting mass) slugs
F_net = F_t + m_net*-g    # (starting net force) lbm



# Matrices: x, v, a, F, m, rho, t, constants
M_n = np.array([[1], [0], [0], [F_net], [m_net], [rho_STP], [0], [1]])  # state matrix

M_s = np.array([[1, t_step, 0, 0, 0, 0, 0, 0],                          # system matrix
                [0, 1, t_step, 0, 0, 0, 0, 0], 
                [0, 0, 0, 1/(M_n[4, 0]), 0, 0, 0, 0],
                [0, -0.5*Cd*A*M_n[5, 0]*abs(M_n[1, 0]), 0, 0, -g, 0, 0, F_t], 
                [0, 0, 0, 0, 1, 0, 0, -m_dot*t_step], 
                [0, 0, 0, 0, 0, 0, 0, rho_STP],
                [0, 0, 0, 0, 0, 0, 1, t_step],
                [0, 0, 0, 0, 0, 0, 0, 1]])

labels = np.array(["Position (ft)",                                     # labels for plots
                   "Velocity (ft/s)", 
                   "Acceleration (ft/s^2)", 
                   "Force (lbf)", 
                   "Mass (slugs)", 
                   "Density (slugs/ft^3)",
                   "Time (sec)"]) 

M_sol = M_n                                                             # solution matrix



# Solver: M_s*V_n = V_n+1
while M_n[0, 0] > 0:

    if M_n[6, 0] >= t_burn:
        M_s[3, 7], M_s[4, 7] = 0, 0

    M_s[3, 1] = -0.5*Cd*A*M_n[5, 0]*abs(M_n[1, 0])
    M_s[5, 7] = (5.1483*10**(-3))/(1 + np.exp((4.74359*10**(-5))*M_n[0, 0] + 0.168201))

    M_n = np.dot(M_s, M_n)

    M_sol = np.hstack([M_sol, M_n])



# Plot all data
for i in range(np.size(M_n) - 2):

    plt.figure()
    plt.plot(M_sol[-2, :], M_sol[i, :])
    plt.xlabel(labels[-1])
    plt.ylabel(labels[i])
    plt.show(block=False)

plt.show()