"""
1D Rocket Flight Simulator
==========================

Author: Ivan Krat 
Initial Release: 2021  
Latest Update: 2025-05-22  
Version: 4.0

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
- Matplotlib: for plotting sim data
- Pandas: for data display

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
v4.0 - Improved drag models and added more features for displaying and exporting simulation data

TODO:
-----
- [x] add drag model
- [x] convert to matrix solver from recursive based solver
- [x] output maximums and minums
- [x] out put csv of all data
- [ ] add multi run and comparison capabiltiy
- [ ] incorperate density and force update back into matrix solver
- [ ] add thrust curves
- [ ] height dependent gravity
- [ ] add parachut staging
- [ ] add super sonic drag modeling
- [x] add skin drag modeling
- [ ] convert to full 3D sim
"""



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



# Constants:
rho_STP = 23.67*(10**-4)  # (desnity at standard tempurature and pressure) slug/ft^3
g = 32.17                 # (gravity) ft/s^2

# Input Variables
D = 11                    # (diameter) in
L = 18                    # (rocket length) ft
Cd = 0.75                 # (ballistic coeff)
Cf = 0.004                # (friction coeff)
F_t = 3000                # (thrust) lbm
m_0 = 225/g               # (dry mass) slugs
m_p = 300/g               # (propellent mass) slugs
m_dot = 10/g              # (mass flow) slugs
t_step = 0.1              # (delta t) sec

# Calculated Values
A = D**2/4*np.pi/144      # (frontal area) ft^2
S = L*D*np.pi/144         # (surface area) ft^2
t_burn = m_p/m_dot        # (burn time) sec
m_net = m_0 + m_p         # (starting mass) slugs
F_net = F_t + m_net*-g    # (starting net force) lbm



# Matrices: x, v, a, F, m, rho, t, constants
M_n = np.array([[1], [0], [0], [F_net], [m_net], [rho_STP], [0], [1]])  # state matrix

M_s = np.array([[1, t_step, 0, 0, 0, 0, 0, 0],                          # system matrix
                [0, 1, t_step, 0, 0, 0, 0, 0], 
                [0, 0, 0, 1/(M_n[4, 0]), 0, 0, 0, 0],
                [0, -0.5*M_n[5, 0]*abs(M_n[1, 0])*(Cd*A + Cf*S), 0, 0, -g, 0, 0, F_t], 
                [0, 0, 0, 0, 1, 0, 0, -m_dot*t_step], 
                [0, 0, 0, 0, 0, 0, 0, rho_STP],
                [0, 0, 0, 0, 0, 0, 1, t_step],
                [0, 0, 0, 0, 0, 0, 0, 1]])

M_sol = M_n                                                             # solution matrix



# Solver: M_s*M_n = M_n+1
while M_n[0, 0] > 0:

    # Thrust Cutoff
    if M_n[6, 0] >= t_burn:
        M_s[3, 7], M_s[4, 7] = 0, 0

    # Update system matrix
    M_s[3, 1] = -0.5*M_n[5, 0]*abs(M_n[1, 0])*(Cd*A + Cf*S)    
    M_s[5, 7] = (5.1483*10**(-3))/(1 + np.exp((4.74359*10**(-5))*M_n[0, 0] + 0.168201))

    # Calculate solution
    M_n = np.dot(M_s, M_n)
    M_sol = np.hstack([M_sol, M_n])



# Display and Save Maximum and Minimum Points
M_sol[4, :] =  M_sol[4, :] * g                    # convert slugs to lbm

max_mins = {                                      # collect all max and min points
    "Data": ["max", "t max", "min", "t min"],
    "Position (ft)": [M_sol[0, :].max(), M_sol[6, M_sol[0, :].argmax()], M_sol[0, :].min(), M_sol[6, M_sol[0, :].argmin()]],
    "Velocity (ft/s)": [M_sol[1, :].max(), M_sol[6, M_sol[1, :].argmax()], M_sol[1, :].min(), M_sol[6, M_sol[1, :].argmin()]],
    "Acceleration (ft/s^2)": [M_sol[2, :].max(), M_sol[6, M_sol[2, :].argmax()], M_sol[2, :].min(), M_sol[6, M_sol[2, :].argmin()]],
    "Force (lbf)": [M_sol[3, :].max(), M_sol[6, M_sol[3, :].argmax()], M_sol[3, :].min(), M_sol[6, M_sol[4, :].argmin()]],
    "Mass (lbm)": [M_sol[4, :].max(), M_sol[6, M_sol[4, :].argmax()], M_sol[4, :].min(), M_sol[6, M_sol[5, :].argmin()]],
    "Density (slugs/ft^3)": [M_sol[5, :].max(), M_sol[6, M_sol[5, :].argmax()], M_sol[5, :].min(), M_sol[6, M_sol[5, :].argmin()]]
}

mm_table = pd.DataFrame(max_mins)
mm_table.to_csv("rocket_maxmin_data.csv", index=False)
print(mm_table) 



# Display and Save Simulation Data
flight_data = {                                   # collect all data for csv
    "Time (sec)": M_sol[6, :],
    "Position (ft)": M_sol[0, :],
    "Velocity (ft/s)": M_sol[1, :],
    "Acceleration (ft/s^2)": M_sol[2, :],
    "Force (lbf)": M_sol[3, :],
    "Mass (lbm)": M_sol[4, :],
    "Density (slugs/ft^3)": M_sol[5, :],
}

fd_table = pd.DataFrame(flight_data)
fd_table.to_csv("rocket_flight_data.csv", index=False)

labels = np.array(["Position (ft)",              # all lables for plots                   
                   "Velocity (ft/s)", 
                   "Acceleration (ft/s^2)", 
                   "Force (lbf)", 
                   "Mass (lbm)", 
                   "Density (slugs/ft^3)",
                   "Time (sec)"])

for i in range(np.size(M_n) - 2):                # plot all variables over time

    plt.figure()
    plt.plot(M_sol[-2, :], M_sol[i, :])
    plt.xlabel(labels[-1])
    plt.ylabel(labels[i])
    plt.show(block=False)

plt.show()