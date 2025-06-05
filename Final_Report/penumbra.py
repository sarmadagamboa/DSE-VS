import numpy as np
import scipy.integrate as integrate

Rm = 3389.500 #km
Rs = 696340.000 #km
dms = 2.28e8 #km
h = 212.478 #km

i = 92.4432 #deg
delta_mars = 25.2 #deg
i_star = np.deg2rad(i) + np.deg2rad(delta_mars) - np.pi/2


### 2D geometry
d_pfoc = Rm/(Rm+Rs) * dms
beta = 2 * np.arcsin(1/d_pfoc * Rm)
delta = np.pi/2 - beta/2
phi = np.arcsin(1/(Rm+h) * Rm)
gamma = np.pi/2 - phi
alpha = 2 * (np.pi - delta - gamma)

epsilon = 2 * np.arctan(Rs / (dms + dms/(Rs-Rm) * Rm))
theta = np.arcsin(Rm/(Rm + h)) - epsilon/2

### 3D geometry
beta_star = np.arcsin(np.sin(np.pi/2-i_star)/np.sin(alpha/2))
a = np.arcsin(np.tan(np.pi/2 - i_star) * (1/np.tan(beta_star)))

percent_in_penumbra = a / np.pi * 100
distance_dip_in_penumbra = (alpha/2 + i_star - np.pi/2) / (alpha/2 - theta/2) * 100
print(f"The percentage of an orbit spent in the penumbra in the worst case is {percent_in_penumbra} %.")
print(f"The furthest dip into the penumbra is {distance_dip_in_penumbra} %, {np.rad2deg(alpha/2 + i_star - np.pi/2)}.")
print(f"a: {np.rad2deg(a)} deg")
print(f"alpha: {np.rad2deg(alpha)} deg")
print(f"theta: {np.rad2deg(theta)} deg")

###integration
#a_star = np.arange(a, 0, 0.001)

def solar_coverage(a_star):
    #c_star = np.arccos(np.multiply(np.cos(a_star), np.cos((np.pi/2) * i_star) * np.ones_like(a_star)))
    c_star = np.arccos(np.cos(a_star) * np.cos((np.pi/2) - i_star))
    #print("------------------")
    #print(f"a_star: {np.rad2deg(a_star)} deg")
    #print(f"c_star: {np.rad2deg(c_star)} deg")
    comp_c_star = alpha/2 - c_star
    #print(f"comp_c_star: {np.rad2deg(comp_c_star)} deg")
    
    cosval = 1 - 4 * comp_c_star / (alpha - theta)
    area_of_sun_covered = (2 * np.arccos(cosval) - np.sin(2 * np.arccos(cosval)))/ (2 * np.pi)
    #print(f"area_of_sun_covered: {area_of_sun_covered}")
    return area_of_sun_covered


int_shadow_half, _ = integrate.quad(solar_coverage, 0, a)
print(f"integrated shadow area half: {int_shadow_half}")
int_shadow = 2 * int_shadow_half
# The factor of 2 accounts for the symmetry of the penumbra on both sides of the orbit.

eq_percent_tot_eclipse = int_shadow / (2*np.pi) * 100

print(f"integrated shadow area: {int_shadow}")
print(f"equivalent percentage of total eclipse: {eq_percent_tot_eclipse} %")
