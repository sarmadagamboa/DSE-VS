import numpy as np
import matplotlib.pyplot as plt

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


# Solar constant
F_s = 1361/1.382**2

# Speed of light
c = 299792458  # m/s

# Reflectance factor
q = 0.6 # ask Bern about this

# Area of the solar panels
A_sp = 4.5 * 2.9

# Area of the side of the spacecraft
A_side = 4.5*1.45

# Area of the antenna
A_a = 1.8**2 * np.pi  # Assuming a circular antenna with diameter 1.8 m

# Surface area
A_s = A_sp + A_side + A_a

# Angle of incidence
i = np.linspace(-0.001, 0.001, 100) # radians, This is the maximum torque case

# Location of center of solar pressure
c_ps = (A_sp*2.9/2 + A_side*(2.9 + 1.45/2) + A_a*(2.9 + 1.45 + 1.8)) / A_s  # Weighted average position of the center of solar pressure, assuming the antenna is perpendicular to the sun

# Location of center of gravity
c_g = 2.9 + 1.45/2  # Assuming the center of gravity is at the center of the spacecraft side, measured from the bottom of the solar panels

F = F_s/c*A_s*(1+q)*np.cos(i)
Tsp = F*(c_ps - c_g)

print(c_ps, c_g)

plt.plot(i, Tsp)
plt.xlabel('Angle of incidence (radians)')
plt.ylabel('Solar Pressure Torque (Nm)')
plt.title('Solar Pressure Torque vs Angle of Incidence')
plt.grid()
plt.show()