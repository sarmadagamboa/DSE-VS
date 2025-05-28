import numpy as np

# Gravitational parameter of Mars (m^3/s^2)
mu = 4.282837e13

# Radius of Mars (m)
r = 3389500

# Orbit radius (m)
R = r + 212000

# Mass moments of inertia (kg*m^2)
I_z = 
I_y = 

# maximum deviation of the Z axis from local vertical (rad)
theta = 

Tg = 3*mu/(2*R**3)*np.abs(I_z - I_y)*np.sin(2*theta)