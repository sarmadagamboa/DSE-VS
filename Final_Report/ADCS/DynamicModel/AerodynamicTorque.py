import numpy as np
from ModelConstants import *

def aerodynamic_torque(C_d, rho, A, r, c_g, c_pa):  # This is constant over one orbit but varies over a martian year

    mu = 42828.3e9  # Gravitational parameter for Mars in m^3/s^2

    # radius of orbit (m)
    R = r + 3.3895e6  # Mars radius + altitude of orbit in m
    
    # Velocity of the spacecraft relative to the atmosphere (m/s)
    V = np.sqrt(mu/R)
    
    F = 0.5*(rho*C_d*A*V**2)
    
    Ta = F*(c_pa - c_g)
    return np.array([0, -Ta[2], Ta[1]])  # Return torque vector in the x, y, z directions

# print(f"Aerodynamic torque: {aerodynamic_torque(C_d, rho, A_aero, r, c_g, c_pa)} Nm, sign according to flight dynamics convention, c_g ofset is {c_pa-c_g}")
