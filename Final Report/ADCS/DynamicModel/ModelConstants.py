import numpy as np


# Moments of inertia for the spacecraft
Ixx = 5000
Iyy = 5000
Izz = 1000
Ixy = 0
Ixz = 0
Iyz = 0

# Initial angular velocity components in rad/s
om_x = 0
om_y = -2*np.pi/(110*60)  # one full revolution over one orbit
om_z = 0

# Minimum torque bits for the control torque in Nm
min_x_bit = 0.0001
min_y_bit = 0.0001
min_z_bit = 0.0001

min_torque = -100
max_torque = 100


##### Disturbance torque constants #####
# Solar radiation pressure torque parameters
q = 0.6 # preliminary reflectance factor, a more comprehensive analysis is needed as it will vary across the spacecraft
A_sp = 4.5 * 2.9 # Area of the solar panels
A_side = 4.5*1.45 # Area of the side of the spacecraft
A_a = 1.8**2 * np.pi  # Area of the antenna # Assuming a circular antenna with diameter 1.8 m
A_s = A_sp + A_side + A_a# Surface area
c_ps = (A_sp*2.9/2 + A_side*(2.9 + 1.45/2) + A_a*(2.9 + 1.45 + 1.8)) / A_s  # Location of center of solar pressure, Weighted average position of the center of solar pressure, assuming the antenna is perpendicular to the sun
c_g_s = 2.9 + 1.45/2  # Location of center of gravity, Assuming the center of gravity is at the center of the spacecraft side, measured from the bottom of the solar panels

# Aerodynamic torque parameters
C_d = 2.6 # Drag coefficient
rho = 1e-11  # Density of the Martian atmosphere in kg/m^3
r = 212e3  # Distance from the center of Mars to the spacecraft in meters
c_g_a = 2.9 + 1.45/2
c_pa = c_g_a + 0.5

# gravity gradient torque parameters
theta = np.pi/2 # Maximum angle between the spacecraft and the local vertical, assuming the spacecraft is in a circular orbit