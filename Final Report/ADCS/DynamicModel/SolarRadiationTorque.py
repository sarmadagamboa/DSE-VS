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


if __name__ == "__main__":
    q = 0.6 # preliminary reflectance factor, a more comprehensive analysis is needed as it will vary across the spacecraft
    # Area of the solar panels
    A_sp = 4.5 * 2.9

    # Area of the side of the spacecraft
    A_side = 4.5*1.45

    # Area of the antenna
    A_a = 1.8**2 * np.pi  # Assuming a circular antenna with diameter 1.8 m

    # Surface area
    A_s = A_sp + A_side + A_a

    # Location of center of solar pressure
    c_ps = (A_sp*2.9/2 + A_side*(2.9 + 1.45/2) + A_a*(2.9 + 1.45 + 1.8)) / A_s  # Weighted average position of the center of solar pressure, assuming the antenna is perpendicular to the sun

    # Location of center of gravity
    c_g = 2.9 + 1.45/2  # Assuming the center of gravity is at the center of the spacecraft side, measured from the bottom of the solar panels

    # Calculate the maximum torque
    Tsp = maximum_torque(q, A_s, c_ps, c_g)
    print(f"Maximum solar radiation torque: {Tsp} Nm, sign according to flight dynamics convention")