import numpy as np
from scipy.constants import G

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

def compute_time_to_distance(delta_v, distance_earth_probe, target_distance):
    """
    Compute the time it takes for a spacecraft to reach a specific distance during the Mars Transfer Injection phase.
    
    Parameters:
    delta_v (float): The delta-v at launch (in km/s).
    distance_earth_probe (float): The distance from Earth to the spacecraft at the start of the mission (in km).
    target_distance (float): The target distance from Earth to probe at a given point (in km).
    
    Returns:
    time_to_reach (float): The time to reach the specified target distance (in days).
    """
    # Convert distances to meters
    distance_earth_probe_m = distance_earth_probe * 1e3  # in meters
    target_distance_m = target_distance * 1e3  # in meters
    
    # Semi-major axis of the transfer orbit (average of Earth's and Mars' orbits)
    a_e = 1 * AU  # Earth orbital radius in km
    a_m = 1.524 * AU  # Mars orbital radius in km
    a_t = (a_e + a_m) / 2  # Semi-major axis of the transfer orbit (Hohmann transfer)

    # Convert to meters
    a_t_m = a_t * 1e3  # semi-major axis in meters

    # Gravitational parameter (mu) of the Sun
    mu = G * M_sun  # in m^3/s^2

    # Orbital period (Kepler's third law)
    T = 2 *np.pi * np.sqrt(a_t_m**3 / mu)  # in seconds
    T_days = T / (60 * 60 * 24)  # Convert to days

    # Using the vis-viva equation to calculate the velocity at any distance r
    def vis_viva(r):
        return np.sqrt(mu * (2/r - 1/a_t_m))

    # Solve for time to reach the target distance by integrating the velocity
    def time_integrand(r):
        return 1 / vis_viva(r)

    # Integrating the time over the path from Earth (distance_earth_probe) to target_distance
    if target_distance_m > distance_earth_probe_m:
        from scipy.integrate import quad
        time_to_reach, _ = quad(time_integrand, distance_earth_probe_m, target_distance_m)
    else:
        return "Target distance must be greater than the starting distance"

    # Convert time from seconds to days
    time_to_reach_days = time_to_reach / (60 * 60 * 24)
    
    return time_to_reach_days


# Example usage: Compute total transfer time and time to reach a specific distance

delta_v_input = float(input("Enter the delta-v at launch (in km/s): "))  # Delta-v in km/s# Distance from Earth to spacecraft in km

# Total transfer time to Mars using Hohmann transfer
transfer_time = compute_transfer_time()
print(f"Total transfer time from Earth to Mars: {transfer_time:.2f} days")

# # Time to reach a specific distance during the MTI phase
# target_distance_input = float(input("Enter the target distance from Earth to spacecraft (in km): "))
# time_to_reach = compute_time_to_distance(delta_v_input, distance_earth_probe_input, target_distance_input)
# print(f"Time to reach {target_distance_input} km: {time_to_reach:.2f} days")
