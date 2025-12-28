"""
1D Rocket Flight Simulator
==========================

Author: Ivan Krat 
Initial Release: 2021  
Latest Update: 2025-05-22  
Version: 4.1

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
- Drag force is calculated using a fixed drag coefficient.
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
- [x] implement Isp
- [ ] add inputs to output csv and prints
- [ ] add multi run and comparison capabiltiy
- [ ] incorperate density and force update back into matrix solver
- [ ] add thrust curves
- [X] add parachut staging
"""



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sympy import true



# outputs
graphs_on = True
max_min_on = True
data_to_csv = True

# Constants:
rho_STP = 23.67*(10**-4)  # (desnity at standard tempurature and pressure) slug/ft^3
g = 32.17                 # (gravity) ft/s^2

# Input Variables
luanch_altitude = 2400    # (altitude of launch site) ft
main_altiude = 1000       # (altitude of main parchute deploy above luanch altitude) ft
D_v = 8                   # (vehical diameter) in
D_d = 4                   # (droug parachute diameter) ft
D_m = 8                   # (main parachute diameter) ft
Cd_v = 0.75               # (vehical drag coeff)
Cd_p = 1.75               # (parachute drag coeff) 
F_t = 1000                # (thrust) lbm
m_0 = 200/g               # (dry mass) lbm/g = slugs
m_p = 100/g               # (propellent mass) lbm/g = slugs
Isp =  250                # (specific impulse) sec
t_step = 0.1              # (delta t) sec

# Calculated Values
A_v = D_v**2/4*np.pi/144  # (vehical frontal area) ft^2
A_d = D_d**2/4*np.pi      # (droug parachute frontal area) ft^2
A_m = D_m**2/4*np.pi      # (main parachute frontal area) ft^2
A = A_v                   # (inital frontal area) ft^2
Cd = Cd_v                 # (inital drag coeff)
m_dot = F_t/Isp/g         # (proppelent mass flow rate) slugs/sec
t_burn =  m_p/m_dot       # (burn time) sec
m_net = m_0 + m_p         # (starting mass) slugs
F_net = F_t + m_net*-g    # (starting net force) lbm
rho = (5.1483*10**(-3))/(1 + np.exp((4.74359*10**(-5))*luanch_altitude + 0.168201))  # (inital atm desnity) slug/ft^3

# (state matrix) x, v, a, F, m, rho, t, constants
M_n = np.array([[luanch_altitude + 1], [0], [0], [F_net], [m_net], [rho], [0], [1]])

# (system matrix)
M_s = np.array([[1, t_step, 0, 0, 0, 0, 0, 0],          # x + v*dt
                [0, 1, t_step, 0, 0, 0, 0, 0],          # v + a*dt
                [0, 0, 0, 1/(M_n[4, 0]), 0, 0, 0, 0],   # a +  F/m
                [0, -0.5*M_n[5, 0]*abs(M_n[1, 0])*(Cd*A), 0, 0, -g, 0, 0, F_t],   # 1/2*rho*A*Cd*v^2 + F_t -g*m
                [0, 0, 0, 0, 1, 0, 0, -m_dot*t_step],   # m - m_dot*dt
                [0, 0, 0, 0, 0, 0, 0, rho],             # rho (updated in solver)
                [0, 0, 0, 0, 0, 0, 1, t_step],          # t + dt
                [0, 0, 0, 0, 0, 0, 0, 1]])              # 1

# (solution matrix) x, v, a, F, m, rho, t, constants
M_sol = M_n



# Solver: M_s*M_n = M_n+1
while M_n[0, 0] > luanch_altitude:

    # Thrust Cutoff
    if M_n[6, 0] >= t_burn:
        M_s[3, 7], M_s[4, 7] = 0, 0

    # Droug Parachute Deploy
    if M_n[1, 0] < 0:
        A = A_d
        Cd = Cd_p
    
    # Main Parachute Deploy
    if M_n[1, 0] < 0 and M_n[0, 0] < (main_altiude+luanch_altitude):
        A = A_m
        Cd = Cd_p

    # Update system matrix
    M_s[3, 1] = -0.5*M_n[5, 0]*abs(M_n[1, 0])*(Cd*A)    
    M_s[5, 7] = (5.1483*10**(-3))/(1 + np.exp((4.74359*10**(-5))*M_n[0, 0] + 0.168201))

    # Calculate solution
    M_n = np.dot(M_s, M_n)
    M_sol = np.hstack([M_sol, M_n])

M_sol[4, :] =  M_sol[4, :] * g                    # convert slugs to lbm



# Display and Save Maximum and Minimum Points
max_mins = {                                      # collect all max and min points
    "Data": ["max", "t max", "min", "t min"],
    "Position (ft)": [M_sol[0, :].max(), M_sol[6, M_sol[0, :].argmax()], M_sol[0, :].min(), M_sol[6, M_sol[0, :].argmin()]],
    "Velocity (ft/s)": [M_sol[1, :].max(), M_sol[6, M_sol[1, :].argmax()], M_sol[1, :].min(), M_sol[6, M_sol[1, :].argmin()]],
    "Acceleration (ft/s^2)": [M_sol[2, :].max(), M_sol[6, M_sol[2, :].argmax()], M_sol[2, :].min(), M_sol[6, M_sol[2, :].argmin()]],
    "Force (lbf)": [M_sol[3, :].max(), M_sol[6, M_sol[3, :].argmax()], M_sol[3, :].min(), M_sol[6, M_sol[4, :].argmin()]],
    "Mass (lbm)": [M_sol[4, :].max(), M_sol[6, M_sol[4, :].argmax()], M_sol[4, :].min(), M_sol[6, M_sol[5, :].argmin()]],
    "Density (slugs/ft^3)": [M_sol[5, :].max(), M_sol[6, M_sol[5, :].argmax()], M_sol[5, :].min(), M_sol[6, M_sol[5, :].argmin()]]
}

if(max_min_on): 
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

if(data_to_csv):
    fd_table = pd.DataFrame(flight_data)
    fd_table.to_csv("rocket_flight_data.csv", index=False)

labels = np.array(["Position (ft)",              # all lables for plots                   
                   "Velocity (ft/s)", 
                   "Acceleration (ft/s^2)", 
                   "Force (lbf)", 
                   "Mass (lbm)", 
                   "Density (slugs/ft^3)",
                   "Time (sec)"])

if(graphs_on):
    for i in range(np.size(M_n) - 2):                # plot all variables over time

        plt.figure()
        plt.plot(M_sol[-2, :], M_sol[i, :])
        plt.xlabel(labels[-1])
        plt.ylabel(labels[i])
        plt.show(block=False)

    plt.show()