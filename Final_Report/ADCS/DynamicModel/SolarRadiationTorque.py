import numpy as np
import matplotlib.pyplot as plt
from ModelConstants import *

## Result of this analysis is that the maximum torque occurs when the sun is perpendicular to the spacecraft

# # Angle of incidence
# i = np.linspace(0, np.pi/2, 100)  # radians

# # Surface area
# A_sc = 4.5*4.35*np.cos(i)

# A_a = 1.8**2*np.pi

# A_s = A_sc + A_a

# x_sc = 4.35/2*np.cos(i)

# x_a = 4.35*np.cos(i) + 1.8

# Delta_x = (A_sc*x_sc + A_a*x_a)/A_s

def maximum_torque(q, A_s, c_ps, c_g):  # This is constant over one orbit but varies over a martian year
    F_s = 1361 / 1.382**2  # Solar constant at Mars
    c = 299792458  # Speed of light in m/s
    i = 0  # Angle of incidence for maximum torque
    F = F_s / c * A_s * (1 + q) * np.cos(i)
    Tsp = F * (c_ps - c_g)
    return np.array([Tsp,0,0])  # Return torque vector in the x, y, z directions
print(f"Maximum solar radiation torque: {maximum_torque(q, A_s, c_ps, c_g_s)} Nm, sign according to flight dynamics convention")