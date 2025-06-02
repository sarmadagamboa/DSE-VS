import numpy as np


def gravity_torque(r, theta, I_y, I_z): # This is constant over one orbit but varies over a martian year
    mu = 42828.3e9  # Gravitational parameter for Mars in m^3/s^2
    R = r + 3.3895e6  # Mars radius + altitude of orbit in m
    Tg = 3*mu/(2*R**3)*np.abs(I_z - I_y)*np.sin(2*theta)
    return np.array([Tg, 0, 0])  # Return torque vector in the x, y, z directions
