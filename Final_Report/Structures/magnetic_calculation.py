import math 

rho_mag_ni = 7180
vol_mag = 0.007 * (110 * 215 + 143.5 * 215)*2*10**(-6)
magnetic_mass = 2 * rho_mag_ni * vol_mag #2 for both spacecraft 

print(magnetic_mass)