import numpy as np

MARS_MU = 4.282837e13  # m^3/s^2
MARS_RADIUS = 3389.5e3  # m

def increase_alt(initial_altitude_m, delta_v1_mps):
    r_pericenter = MARS_RADIUS + initial_altitude_m
    v_circular = np.sqrt(MARS_MU / r_pericenter)

    # After first burn at periapsis
    v_pericenter_new = v_circular + delta_v1_mps

    # New semi-major axis
    a_transfer = 1 / (2 / r_pericenter - v_pericenter_new**2 / MARS_MU)

    # Apoapsis of the new elliptical orbit
    r_apocenter = 2 * a_transfer - r_pericenter
    apocenter_altitude = r_apocenter - MARS_RADIUS

    # Velocity at apoapsis of elliptical orbit
    v_apocenter_elliptical = np.sqrt(MARS_MU * (2 / r_apocenter - 1 / a_transfer))

    # Velocity needed to circularize at apoapsis
    v_circular_apocenter = np.sqrt(MARS_MU / r_apocenter)

    delta_v2 = v_circular_apocenter - v_apocenter_elliptical

    return {
        "delta_v2_mps": delta_v2,
        "final_circular_orbit_altitude_m": apocenter_altitude
    }


def delta_v_alt(initial_altitude_m, final_altitude_m):
    r1 = MARS_RADIUS + initial_altitude_m
    r2 = MARS_RADIUS + final_altitude_m

    # Circular velocities at each orbit
    v1 = np.sqrt(MARS_MU / r1)
    v2 = np.sqrt(MARS_MU / r2)

    # Semi-major axis of transfer orbit
    a_transfer = (r1 + r2) / 2

    # Velocities on transfer orbit at r1 and r2
    v_transfer1 = np.sqrt(MARS_MU * (2 / r1 - 1 / a_transfer))
    v_transfer2 = np.sqrt(MARS_MU * (2 / r2 - 1 / a_transfer))

    # Delta-vs for each burn
    delta_v1 = v_transfer1 - v1
    delta_v2 = v2 - v_transfer2

    total_delta_v = np.abs(delta_v1) + np.abs(delta_v2)

    return {
        "delta_v1_mps": delta_v1,
        "delta_v2_mps": delta_v2,
        "total_delta_v_mps": total_delta_v
    }


# altitude = 212e3  
# delta_v1 = 0.004 #m/s
# result = increase_alt(altitude, delta_v1)  # 400 km start, 100 m/s first burn

# print(f"Second Burn Δv: {result['delta_v2_mps']:.2f} m/s")
# print(f"Final Circular Orbit Altitude: {result['final_circular_orbit_altitude_m']/1000:.2f} km")


altitude1 = 212.5e3-0.313e3
altitude2 = 212.5e3+0.313e3  
result = delta_v_alt(altitude1, altitude2)
print(f"Δv1 (Transfer Burn): {result['delta_v1_mps']:.2f} m/s")
print(f"Δv2 (Circularization): {result['delta_v2_mps']:.2f} m/s")
print(f"Total Δv: {result['total_delta_v_mps']:.2f} m/s")
