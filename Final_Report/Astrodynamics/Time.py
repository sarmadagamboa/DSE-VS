import numpy as np
from scipy.constants import G
from marsdensity import plot_mars_density
# Constants
M_sun = 1.989e30  # Mass of the Sun in kg
AU = 1.496e8  # Astronomical Unit in km

def compute_transfer_time():
    """
    Compute the time it takes for a spacecraft to cover a certain distance during the Mars Transfer Injection phase.

    Parameters:
    delta_v (float): The delta-v at launch (in km/s).
    distance_earth_probe (float): The distance from Earth to the spacecraft at the start of the mission (in km).

    Returns:
    time_to_transfer (float): The time to transfer from Earth to Mars using Hohmann transfer (in days).
    """
    # Convert distance to meters
    
    # Semi-major axis of the transfer orbit (average of Earth's and Mars' orbits)
    a_e = 1 * AU  # Earth orbital radius in km
    a_m = 1.524 * AU  # Mars orbital radius in km
    a_t = (a_e + a_m) / 2  # Semi-major axis of the transfer orbit (Hohmann transfer)

    # Convert to meters
    a_t_m = a_t * 1e3  # semi-major axis in meters

    # Gravitational parameter (mu) of the Sun
    mu = G * M_sun  # in m^3/s^2

    # Period of the transfer orbit (Kepler's third law)
    T = 2 * np.pi * np.sqrt(a_t_m**3 / mu)  # in seconds
    T_days = T / (60 * 60 * 24)  # Convert to days

    # Using the spacecraft velocity at launch is the orbital velocity at Earth's orbit # Convert to meters

    # Using the two-impulse Hohmann transfer delta-v to compute the total time
    # The time to transfer from Earth to Mars is roughly half the orbital period
    time_to_transfer = T_days / 2  # As it's a half orbit transfer

    return time_to_transfer

def orbital_decay(alt,start_apocenter = 45000, year= 2, spacecraft_mass=595*2, drag_coefficient=2.2, cross_sectional_area=1.54):
    """
    Calculate the orbital decay due to atmospheric drag over a specified time period.
    """
    H = 11e3
    R_mars = 3390e3
    mu = 4.282837e13  # Mars gravitational parameter in m^3/s^2

    results = plot_mars_density(year)
    heights = results['height']
    densities = results['max_density']
    
    idx = np.abs(heights - alt).argmin()
    closest_alt = heights[idx]
    rho = densities[idx]
    print(rho)
    if np.isnan(rho):
        raise ValueError(f"No atmospheric density data available for altitude {alt} km.")

    r = R_mars +closest_alt*1000
    apocenter = R_mars + start_apocenter * 1000  # Convert to meters
    a_current = (r + apocenter)/2 # Semi-major axis in meters
    v_pericenter = np.sqrt(mu * (2 / r - 1 / a_current))
    drag_force = 0.5 * rho * v_pericenter**2 * drag_coefficient * cross_sectional_area
    print(drag_force)

    s = np.pi * r  # Circumference of the orbit
    delta_E = -drag_force * s 
    E = -spacecraft_mass * mu / (2*a_current)  # Initial orbital energy
    E_final = E + delta_E
    a_final = -spacecraft_mass * mu / (2 * E_final)
    #print(r_final)  # New radius after decay
    delta_r = (-a_current+a_final)/1000
    return delta_r

lifetime = 1.88  # Operational lifetime in years
alt = 100
# Total transfer time to Mars using Hohmann transfer
transfer_time = compute_transfer_time()
print(f"Total transfer time from Earth to Mars: {transfer_time:.2f} days")
total_decay = orbital_decay(
    alt=alt,
    year=2,  # 2 years of decay
    spacecraft_mass=595*2,  # kg
    drag_coefficient=2.2,
    cross_sectional_area=1.54,  # m/s orbital velocity   # years of operation before decay
)
print(f"Total orbital decay at {alt} km altitude: {total_decay:.2f} km")
