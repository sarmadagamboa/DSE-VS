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

def calculate_aerobraking_time_np(initial_apocenter_altitude_km=212.5, spacecraft_mass_kg=654.731,drag_coefficient=2.6, cross_sectional_area_m2=1.45):
    print("--- Mars Aerobraking Time Calculator (using NumPy) ---")
    print("This program estimates the aerobraking duration on Mars.")
    print("It uses a simplified exponential atmospheric model and assumes a constant pericenter.")
    print("Note: Aerobraking at 212.5 km pericenter is extremely slow due to very thin atmosphere.")
    print("Actual aerobraking maneuvers typically occur at much lower altitudes (e.g., 100-150 km).")
    print("-" * 40)

    G = 6.67430e-11
    M_mars = 6.4169e23
    MU_MARS = G * M_mars
    R_MARS = 3389.5e3

    PERICENTER_DENSITY = 1.797e-11
    H_SCALE = 11.1e3


    PERICENTER_ALTITUDE_KM = 212.5
    TARGET_APOCENTER_ALTITUDE_KM = 212.5

    
    initial_apocenter_radius = R_MARS + initial_apocenter_altitude_km * 1000
    pericenter_radius = R_MARS + PERICENTER_ALTITUDE_KM * 1000
    target_apocenter_radius = R_MARS + TARGET_APOCENTER_ALTITUDE_KM * 1000

    current_apocenter_altitude = initial_apocenter_altitude_km
    total_time_seconds = 0.0
    num_orbits = 0

    print("\nStarting aerobraking simulation...")
    MAX_ORBITS = 10_000_000
    APOCENTER_TOLERANCE_KM = 0.1

    while current_apocenter_altitude > (TARGET_APOCENTER_ALTITUDE_KM + APOCENTER_TOLERANCE_KM) and num_orbits < MAX_ORBITS:
        current_apocenter_radius = R_MARS + current_apocenter_altitude * 1000
        
        current_semi_major_axis = (pericenter_radius + current_apocenter_radius) / 2

        if current_semi_major_axis <= 0:
            print("Error: Semi-major axis became non-positive. Aborting simulation.")
            break

        orbital_period_seconds = 2 * np.pi * np.sqrt(current_semi_major_axis**3 / MU_MARS)

        velocity_at_pericenter = np.sqrt(MU_MARS * (2 / pericenter_radius - 1 / current_semi_major_axis))

        density_at_pericenter = PERICENTER_DENSITY

        delta_semi_major_axis_per_orbit = - (drag_coefficient * cross_sectional_area_m2 / spacecraft_mass_kg) * \
                                            density_at_pericenter * velocity_at_pericenter**2 * \
                                            H_SCALE * (2 * current_semi_major_axis**2 / MU_MARS)
        
        delta_apocenter_radius_per_pass = 2 * delta_semi_major_axis_per_orbit

        delta_apocenter_altitude_per_pass_km = delta_apocenter_radius_per_pass / 1000
        print(f"THIS IS THE DELTA ALTITUDE:{delta_apocenter_altitude_per_pass_km} km")
        current_apocenter_altitude += delta_apocenter_altitude_per_pass_km

        total_time_seconds += orbital_period_seconds
        num_orbits += 1

        if num_orbits % 10000 == 0 or (num_orbits < 1000 and num_orbits % 100 == 0) :
            print(f"  Orbit: {num_orbits:,.0f}, Current Apocenter: {current_apocenter_altitude:.2f} km, Time: {total_time_seconds / (24 * 3600):.2f} days")

    total_time_days = total_time_seconds / (24 * 3600)
    total_time_years = total_time_days / 365.25

    print("\n--- Aerobraking Simulation Complete ---")
    if num_orbits >= MAX_ORBITS:
        print(f"Simulation stopped after {MAX_ORBITS:,.0f} orbits (Max limit reached).")
        print("This indicates a very long aerobraking duration or insufficient drag for the target.")
    else:
        print("Target apocenter altitude reached.")

    print(f"\nInitial Apocenter Altitude: {initial_apocenter_altitude_km:.2f} km")
    print(f"Pericenter Altitude: {PERICENTER_ALTITUDE_KM:.2f} km (Constant)")
    print(f"Final Apocenter Altitude: {current_apocenter_altitude:.2f} km")
    print(f"Spacecraft Mass: {spacecraft_mass_kg:.2f} kg")
    print(f"Cross-sectional Area: {cross_sectional_area_m2:.2f} m^2")
    print(f"Drag Coefficient: {drag_coefficient:.2f}")

    print(f"\nTotal Number of Orbits: {num_orbits:,.0f}")
    print(f"Total Aerobraking Time: {total_time_days:,.2f} days")
    print(f"                          ({total_time_years:,.2f} years)")

calculate_aerobraking_time_np()

# Total transfer time to Mars using Hohmann transfer
transfer_time = compute_transfer_time()
print(f"Total transfer time from Earth to Mars: {transfer_time:.2f} days")

