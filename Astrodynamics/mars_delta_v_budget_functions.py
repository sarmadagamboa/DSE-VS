
import numpy as np
from astropy.constants import G, M_earth, R_earth
from astropy import units as u

# Constants
M_mars = 6.4171e23 * u.kg
R_mars = 3389.5 * u.km
d_MS = 228e6  # km
d_ES = 149e6  # km
mu_earth = (G * M_earth).to(u.km**3 / u.s**2).value
mu_mars = (G * M_mars).to(u.km**3 / u.s**2).value
r_earth = R_earth.to(u.km).value
r_mars = R_mars.to(u.km).value
mu_sun = 1.32712440018e11  # km^3/s^2

# Mission parameters
leo_alt = 300
mars_orbit_alt = 300
delta_v_per_year = 60  # m/s
leo_radius = r_earth + leo_alt
mars_orbit_radius = r_mars + mars_orbit_alt
v_leo = np.sqrt(mu_earth / leo_radius)
v_mars_circ = np.sqrt(mu_mars / mars_orbit_radius)


def compute_launch_to_leo():
    return np.sqrt(mu_earth / leo_radius)

def compute_transfer_injection():
    transfer_a = (leo_radius + (d_MS + d_ES) + mars_orbit_radius) / 2
    v_earth_orbit = np.sqrt(mu_sun / d_ES)
    v_hel_LEO = np.sqrt(mu_sun * (2 / d_ES - 1 / transfer_a))
    v_inf_leo = v_hel_LEO - v_earth_orbit
    v_per_LEO = np.sqrt((2 * mu_earth / leo_radius) + v_inf_leo**2)
    delta_v = np.abs(v_per_LEO - v_leo)
    eccentricity = 1 - ((leo_radius + d_ES) / transfer_a)
    return delta_v, transfer_a, eccentricity, v_inf_leo


def compute_mars_inclination_change(v_orbit):

    delta_v_incl = np.sqrt(2 * v_orbit**2 * (1 - np.cos(np.radians(93))))
    return delta_v_incl, r_apo_mars

def compute_mars_circularization(v_inf_mars):
    v_apo_mars = np.sqrt((2 * mu_mars / mars_orbit_radius) + v_inf_mars**2)
    return np.abs(v_mars_circ - v_apo_mars)

def compute_station_keeping(years=4.5):
    return delta_v_per_year * years / 1000  # km/s

def compute_deorbit():
    return 60 / 1000  # km/s

# Run all calculations
e_mars = 0.8537
a_mars = mars_orbit_radius / (1 - e_mars)
r_apo_mars = 2 * a_mars - mars_orbit_radius
v_mars_apo = np.sqrt(mu_mars / r_apo_mars)



delta_v_launch = compute_launch_to_leo()
delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection()
v_hel_mars = np.sqrt(mu_sun * (2 / d_MS - 1 / transfer_a))
v_mars_orbit = np.sqrt(mu_sun / d_MS)
v_inf_mars = v_hel_mars - v_mars_orbit

# r_apo_mars = compute_mars_inclination_change(v_inf_mars)

# delta_v_2 = 0  # Aerobraking assumed
use_aerobraking = True  

if use_aerobraking:
    delta_v_2 = 0
    delta_v_inclination = compute_mars_inclination_change(v_mars_apo)[0]
    print("Aerobraking enabled — circularization ΔV saved.")
else:
    delta_v_2 = compute_mars_circularization(v_inf_mars)
    delta_v_inclination = compute_mars_inclination_change(v_mars_circ)[0]
    print("Aerobraking disabled — full circularization burn required.")



delta_v_station_keeping = compute_station_keeping()
delta_v_deorbit = compute_deorbit()

# Set to False if no aerobraking

# Total ΔV
total_delta_v = (
    delta_v_launch +
    delta_v_1 +
    delta_v_inclination +
    delta_v_2 +
    delta_v_station_keeping +
    delta_v_deorbit
)

# Print Results
print(f'Eccentricity of transfer orbit: {e:.4f}')
print("=== Mars Mission ΔV Budget ===")
print(f'ΔV: Launch to LEO              = {delta_v_launch:.3f} km/s')
print(f'ΔV: LEO to Mars Transfer       = {delta_v_1:.3f} km/s')
print(f'ΔV: Inclination Change         = {delta_v_inclination:.3f} km/s')
print(f'ΔV: MARS Circularization       = {delta_v_2:.3f} km/s')
print(f'ΔV: Station Keeping (4.5 yrs)  = {delta_v_station_keeping:.3f} km/s')
print(f'ΔV: End-of-Life Deorbit        = {delta_v_deorbit:.3f} km/s')
print(f'------------------------------------------')
print(f'Total Mission ΔV               = {total_delta_v:.3f} km/s')
