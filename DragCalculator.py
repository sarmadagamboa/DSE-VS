import numpy as np
import matplotlib.pyplot as plt

# Density at 175 km altitude https://ui.adsabs.harvard.edu/abs/2001JGR...10623349T/abstract
rho = 1.8e-11 # kg/m^3 

GM = 0.042828e6 # km^3/s^2
R = 3.3895e6 # m
h = 250e3 # m
r = R + h # m
r = r / 1000 # km
A = 12.5 + 2.4*2.4
C_d = 2.6

v = np.sqrt(GM/r) # km/s


D = 0.5*rho*(v*1000)**2 * A * C_d # N

lifetime = 4.5 # years

m = 600 # kg

deltaV = D/m*lifetime*365*24*3600 # m/s

print(f'for a spacecraft at {h/1000} km, traveling at {v} km/s, the drag force is {D} N. Over a lifetime of {lifetime} years, the drag amounts to a ΔV of {deltaV} m/s.')