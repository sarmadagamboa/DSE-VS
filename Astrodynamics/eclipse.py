import numpy as np
from datetime import datetime, timedelta
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.astro.time_conversion import calendar_date_to_julian_day, julian_day_to_calendar_date
from tudatpy.kernel.astro.epoch import epoch_from_julian_day

# Load SPICE kernels
spice.load_standard_kernels()

# ---------------------
# Simulation Settings
# ---------------------
start_utc = datetime.utcnow()
duration_minutes = 120
step_seconds = 30

# Mars central body and relevant celestial objects
central_body = 'Mars'
bodies_to_create = ['Mars', 'Sun']
global_frame_origin = 'Mars'
global_frame_orientation = 'J2000'

# Create system of bodies
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create,
    global_frame_origin,
    global_frame_orientation
)
bodies = environment_setup.create_system_of_bodies(body_settings)

# Mars gravitational parameter
mu_mars = bodies.get_body('Mars').gravity_field_model.get_gravitational_parameter()

# Mars SSO-like orbital elements (approx. 400 km, 93° inclination)
mars_radius = 3396.2e3  # meters
altitude = 400e3
semi_major_axis = mars_radius + altitude
eccentricity = 0.0
inclination = np.deg2rad(93.0)
RAAN = np.deg2rad(0.0)
arg_periapsis = 0.0
true_anomaly = 0.0

kepler_elements = [semi_major_axis, eccentricity, inclination, RAAN, arg_periapsis, true_anomaly]
initial_state = element_conversion.keplerian_to_cartesian(kepler_elements, mu_mars)

# Time vector in Julian Days
times = [calendar_date_to_julian_day(start_utc + timedelta(seconds=i))
         for i in range(0, duration_minutes * 60, step_seconds)]

# ---------------------
# Eclipse Geometry
# ---------------------
sun_radius = 696340e3  # meters

def compute_eclipse_type(spacecraft_pos, sun_pos, planet_radius, sun_radius):
    """Returns UMBRA, PENUMBRA, or FULL using 3D shadow cone logic"""
    sun_to_planet = -sun_pos
    sun_to_sat = spacecraft_pos - sun_pos

    d_sp = np.linalg.norm(sun_to_planet)
    d_ss = np.linalg.norm(sun_to_sat)

    alpha = np.arcsin(planet_radius / d_sp)  # Angular radius of planet
    beta = np.arcsin(sun_radius / d_ss)      # Angular radius of sun
    angle = np.arccos(np.dot(sun_to_planet, sun_to_sat) / (d_sp * d_ss))

    if angle < (alpha - beta):
        return "UMBRA"
    elif angle < (alpha + beta):
        return "PENUMBRA"
    else:
        return "FULL"

# ---------------------
# Run Eclipse Check
# ---------------------
eclipse_log = []
previous_status = "FULL"

for jd in times:
    t = epoch_from_julian_day(jd)

    # Get Mars-centered inertial position of the Sun
    sun_pos = spice.get_body_cartesian_state_at_epoch(
        target_body_name='Sun',
        observer_body_name='Mars',
        reference_frame_name=global_frame_orientation,
        aberration_corrections='none',
        ephemeris_time=t
    )[:3]

    # Spacecraft position (inertial, constant orbital elements for now)
    spacecraft_pos = initial_state[:3]

    eclipse_type = compute_eclipse_type(
        spacecraft_pos,
        sun_pos,
        planet_radius=mars_radius,
        sun_radius=sun_radius
    )

    # Detect transitions
    if eclipse_type != previous_status:
        eclipse_log.append((julian_day_to_calendar_date(jd), f"ENTER {eclipse_type}"))
        previous_status = eclipse_type

# ---------------------
# Display Results
# ---------------------
print(f"Eclipse events for Mars orbiter from {start_utc} UTC:\n")
for time, status in eclipse_log:
    print(f"{time}: {status}")
