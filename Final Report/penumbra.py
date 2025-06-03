import math as m
import numpy as np

Rm = 10 #update
Rs =1000 #update
dms =1000000 #update
h = 1 #update

i = 3 #update
delta_mars = 25 #update
i_star = i + delta_mars


### 2D geometry
d_pfoc = Rm/(Rm+Rs) * dms
beta = 2 * m.asin(1/d_pfoc * Rm)
delta = m.pi/2 - beta/2
phi = m.asin(1/(Rm+h) * Rm)
gamma = m.pi/2 - phi
alpha = 2 * (m.pi - delta - gamma)


### 3D geometry
beta_star = m.asin(m.sin(m.pi/2-i_star)/m.sin(alpha/2))
a = m.asin(m.tan(m.pi/2 - i_star) * (1/m.tan(beta_star)))

percent_in_penumbra = a / m.pi * 100
distance_dip_in_penumbra = alpha/2 + i_star - m.pi/2
print(f"The percentage of an orbit spent in the penumbra in the worst case is {percent_in_penumbra}%.")
print(f"The furthest dip into the penumbra is {distance_dip_in_penumbra}%.")
print(f"a: {a}")

###integration
a_star = np.arange(a, 0, 0.001)

c_star = np.arccos(np.multiply(np.cos(a_star), np.cos((np.pi/2) * i_star) * np.ones_like(a_star)))
area_of_sun_covered = (2 * np.arccos(1-c_star/alpha) - )