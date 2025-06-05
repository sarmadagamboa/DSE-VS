import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.constants import R_sun
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from poliastro.bodies import Mars
from poliastro.twobody import Orbit
from poliastro.frames import Planes

def rotate_about_x(vec, angle_rad):
    """Rotate vector about X-axis by angle_rad."""
    R = np.array([
        [1, 0, 0],
        [0, np.cos(angle_rad), -np.sin(angle_rad)],
        [0, np.sin(angle_rad),  np.cos(angle_rad)]
    ])
    return R @ vec

def get_solstice_vectors(altitude_km, epoch_str="2025-07-12 12:00:00"):
    epoch = Time(epoch_str, scale="tdb")

    # Get Mars and Sun positions
    with solar_system_ephemeris.set("jpl"):
        r_mars = get_body_barycentric("mars", epoch).xyz.to(u.km).value
        r_sun = get_body_barycentric("sun", epoch).xyz.to(u.km).value

    # Mars obliquity (tilt)
    obliquity_deg = 25.19
    obliquity_rad = np.radians(obliquity_deg)

    # Mars→Sun vector in inertial frame
    r_msun_inertial = r_sun - r_mars
    r_msun_equatorial = rotate_about_x(r_msun_inertial, -obliquity_rad)  # rotate Sun vector into Mars equatorial frame

    # Orbit in Mars equatorial plane, 93° SSO
    orbit = Orbit.circular(
        Mars,
        alt=altitude_km * u.km,
        inc=93 * u.deg,
        epoch=epoch,
        plane=Planes.EARTH_EQUATOR  # use Mars equator as base
    )
    r_sc = orbit.rv()[0].to(u.km).value  # spacecraft position

    r_os = r_sc - r_sun
    r_ms = r_sc - r_mars

    return r_os, r_ms, r_sc, np.linalg.norm(r_msun_equatorial)

def compute_eclipse_duration_conical_model_solstice(altitude_km, epoch_str="2025-07-12 12:00:00"):
    # Constants
    R_M = 3390.0  # Mars radius in km
    R_S = R_sun.to(u.km).value  # Sun radius in km
    mu_mars = 42828.3  # Mars gravitational parameter, km^3/s^2

    r_os, r_ms, r_sc, r_os_mag = get_solstice_vectors(altitude_km, epoch_str)

    r_hat_os = r_os / np.linalg.norm(r_os)
    r_hat_ms = r_ms / np.linalg.norm(r_ms)

    dot_product = np.dot(r_hat_os, r_ms)
    r_ps = dot_product * r_hat_os
    r_ps_mag = np.linalg.norm(r_ps)

    r_mp = r_ms - r_ps
    h_g = np.linalg.norm(r_mp) - R_M

    R_p = (r_ps_mag / r_os_mag) * R_S

    h_c = h_g
    h_t = h_c + R_p
    h_b = h_c - R_p

    eta = h_c / (h_c - h_b)

    if eta < -1:
        eclipse_type = "Full eclipse (umbra)"
    elif -1 <= eta <= 1:
        eclipse_type = "Partial eclipse (penumbra)"
    else:
        eclipse_type = "Fully sunlit"

    r_sat = R_M + altitude_km
    theta_penumbra = np.arcsin(R_M / r_sat) + np.arcsin(R_S / r_os_mag)
    orbit_period = 2 * np.pi * np.sqrt(r_sat**3 / mu_mars)
    eclipse_duration = (2 * theta_penumbra / (2 * np.pi)) * orbit_period

    if abs(eta) <= 1:
        f_g = 1 - (1 / np.pi) * np.arccos(eta) + (eta / np.pi) * np.sqrt(1 - eta**2)
    elif eta < -1:
        f_g = 1.0
    else:
        f_g = 0.0

    return {
        "altitude_km": altitude_km,
        "eta": eta,
        "eclipse_type": eclipse_type,
        "eclipse_fraction": f_g,
        "eclipse_time_s": eclipse_duration if f_g > 0 else 0,
        "orbit_period_min": orbit_period / 60,
        "epoch": epoch_str
    }

# Example usage
result = compute_eclipse_duration_conical_model_solstice(altitude_km=212)

print("--- Mars Eclipse Summary (Solstice) ---")
for k, v in result.items():
    if isinstance(v, float):
        print(f"{k}: {v:.6f}")
    else:
        print(f"{k}: {v}")
