# delta V budget
from astropy import units as u
from astropy.constants import G, M_earth, R_earth
import numpy as np


# Add these lines manually to define Mars constants
M_mars = 6.4171e23 * u.kg     # Mars mass
R_mars = 3389.5 * u.km        # Mars radius (mean)
d_MS = 228e6 
d_ES = 149e6 
# Constants
mu_earth = (G * M_earth).to(u.km**3 / u.s**2).value
mu_mars = (G * M_mars).to(u.km**3 / u.s**2).value
r_earth = R_earth.to(u.km).value
r_mars = R_mars.to(u.km).value
mu_sun = 1.32712440018e11  # km^3/s^2 (Sun)
delta_v_per_year = 60 #m/s
# Mission Parameters
leo_alt = 300  # km
leo_radius = r_earth + leo_alt  # km

mars_orbit_alt = 500  # km
mars_orbit_radius = r_mars + mars_orbit_alt  # km

mars_sma = (d_MS + d_ES) # km (major axis of Mars orbit)
#Phase 0: Launch - LEO
delta_v_launch = np.sqrt(mu_earth/leo_radius)

#Planetary velocities
v_earth_orbit = np.sqrt(mu_sun/d_ES)
v_mars_orbit = np.sqrt(mu_sun / d_MS)

# Earth LEO Orbit - MTI
transfer_a = ((leo_radius + mars_sma + mars_orbit_radius) / 2)  # km
e = 1-((leo_radius+d_ES)/transfer_a)
print(transfer_a)
v_leo = np.sqrt(mu_earth / leo_radius)
v_hel_LEO = np.sqrt(mu_sun * (2 /d_ES - 1 / transfer_a))
v_inf_leo = v_hel_LEO - v_earth_orbit
v_per_LEO = np.sqrt((2*mu_earth/leo_radius) + v_inf_leo**2)
delta_v_1 =  np.abs(v_per_LEO - v_leo)

# Mars hyperbolic velocity
v_mars = np.sqrt(mu_mars/mars_orbit_radius) #this is the velocity at the desired altitude
v_hel_mars = np.sqrt(mu_sun * (2 /d_MS - 1 / transfer_a))
v_inf_mars = v_hel_mars - v_mars_orbit



# Phase 4: Circularization at altitude
v_apo_mars = np.sqrt((2*mu_mars/mars_orbit_radius) + v_inf_leo**2)
delta_v_inclination = np.sqrt(2*v_mars**(2) *(1-np.cos(93)))
delta_v_2= np.abs(v_mars - v_apo_mars)

# Phase 5: Station Keeping (15 m/s per year for 4.5 years)

delta_v_station_keeping = delta_v_per_year * 4.5 / 1000  # km/s

# Phase 6: End-of-Life Decommissioning Burn
delta_v_deorbit = 60/1000  # km/s

# Total ΔV
total_delta_v = (
    delta_v_launch+
    delta_v_1+
    delta_v_inclination+
    delta_v_2+
    delta_v_station_keeping+
    delta_v_deorbit
)

# Print Results
print(f'eccentricity of transfer orbit: {e}')
print("=== Mars Mission ΔV Budget ===")
print(f'ΔV: launch to LEO = {delta_v_launch:.3f} km/s')
print(f'ΔV: from LEO - MARS = {delta_v_1:.3f} km/s')
print(f'ΔV: orbit inclination = {delta_v_inclination:.3f} km/s')
print(f'ΔV: MARS circularize= {delta_v_2:.3f} km/s')
print(f"ΔV: Station Keeping      = {delta_v_station_keeping:.3f} km/s")
print(f"ΔV: End-of-Life Deorbit             = {delta_v_deorbit:.3f} km/s")
print(f"------------------------------------------")
print(f"Total Mission ΔV                    = {total_delta_v:.3f} km/s")
