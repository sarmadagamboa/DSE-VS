from astropy import units as u
from astropy.time import Time, TimeDelta
from poliastro.bodies import Sun, Mars
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import propagate
from poliastro.ephem import build_ephem_interpolant
from poliastro.maneuver import Maneuver
from poliastro.frames import Planes
from poliastro.util import time_range
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import get_body_barycentric





def compute_eclipse_duration():
    # Define simulation time
    start_time = Time("2025-01-01 00:00:00", scale="utc")

    # Mars parameters
    mars_radius = Mars.R.to(u.km).value  # ~3396.2 km

    # Sun-synchronous orbit approximations: dawn-dusk => near-polar, ~low orbit
    # Let's assume ~200 km altitude circular orbit
    alt = 470 * u.km
    a = mars_radius * u.km + alt  # semi-major axis
    ecc = 0.0 * u.one
    inc = 92 * u.deg  # near polar orbit (dawn-dusk)
    raan = 0 * u.deg
    argp = 0 * u.deg
    nu = 0 * u.deg

    # Create orbit
    orb = Orbit.from_classical(Mars, a, ecc, inc, raan, argp, nu, epoch=start_time)

    # Determine period and create time range for one orbit
    period = orb.period
    print(f"Orbit period: {period.to(u.s)/60} minutes")
    spacing = TimeDelta((period / 500).to(u.s))  # Divide the orbit into 500 intervals
    times = time_range(start=start_time, spacing=spacing, periods=500)

    # Get Sun position relative to Mars at each time
    sun_positions = []
    for t in times:
        mars_to_sun = get_body_barycentric("sun", t) - get_body_barycentric("mars", t)
        sun_positions.append(mars_to_sun.xyz.to(u.km).value)


    sun_positions = np.array(sun_positions)  # shape (N, 3)

    # Propagate spacecraft and compute positions in Mars-centered frame
    sc_positions = []
    for t in times:
        sc_pos = orb.propagate(t - orb.epoch).rv()[0]
        sc_positions.append(sc_pos.to(u.km).value)

    sc_positions = np.array(sc_positions)

    # Check for eclipse: simple geometric shadow model (umbra cone)
    eclipse_flags = []
    for sc_pos, sun_pos in zip(sc_positions, sun_positions):
        # Vector from Mars to Sun and to spacecraft
        v_sun = sun_pos
        v_sat = sc_pos

        # Angle between Sun and spacecraft vectors
        angle = np.arccos(np.dot(v_sun, v_sat) / (np.linalg.norm(v_sun) * np.linalg.norm(v_sat)))

        # Shadow cone approximation: Mars angular radius from spacecraft
        mars_angular_radius = np.arcsin(mars_radius / np.linalg.norm(v_sat))

        # Eclipse if spacecraft behind Mars relative to Sun (angle < Mars angular radius)
        eclipse = angle < mars_angular_radius
        eclipse_flags.append(eclipse)

    eclipse_flags = np.array(eclipse_flags)

    # Compute eclipse duration
    dt = (period / len(times)).to(u.s)
    total_eclipse_time = np.sum(eclipse_flags) * dt

    print(f"Total eclipse duration over one orbit: {total_eclipse_time:.2f}")

    # Optional: plot geometry
    plt.figure()
    plt.plot(times.datetime, eclipse_flags.astype(int), drawstyle="steps-post")
    plt.ylabel("In Eclipse")
    plt.xlabel("Time (UTC)")
    plt.title("Eclipse Events Over One Orbit")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    compute_eclipse_duration()
