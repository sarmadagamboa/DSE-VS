import numpy as np

# Drag coefficient
C_d = 2.6  # Should calculate this

# Density of the atmosphere at Mars (kg/m^3)
rho = 

# Surface area of the spacecraft (m^2)
A = 

# Velocity of the spacecraft relative to the atmosphere (m/s)
V =

# Center of aerodynamic pressure (m)
c_pa = 

# Center of gravity (m)
c_g =

F = 0.5*(rho*C_d*A*V**2)
Ta = F*(c_pa - c_g)