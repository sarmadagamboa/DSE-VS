import numpy as np
from astropy.constants import G
from astropy import units as u

def calc_2d_elipse(h): # h in km

    R_mars = 3389.5 # km
    M_mars = 6.4171e23 * u.kg
    mu_mars = (G * M_mars).to(u.km**3 / u.s**2).value

    theta = np.degrees(np.arccos((R_mars/(R_mars + h))))

    eclipse_percent = (180-2*theta) / 360 * 100

    eclipse_time = eclipse_percent * 2 * np.pi * np.sqrt((R_mars + h)**3 / mu_mars)

    return eclipse_time, eclipse_percent


# def mars_density_model(altitude, Ls, latitude, local_time, F10.7, dust_tau):
#     """
#     Returns density (kg/m³), composition (CO₂, O, He, H₂), and scale height.
#     """
#     # Base density from Mars-GRAM 2021
#     rho_base = mars_gram(altitude, Ls, latitude)
    
#     # Thermospheric corrections (MAVEN data)
#     rho_thermo = rho_base * (1 + 0.3 * np.sin(local_time * np.pi/12))
    
#     # Solar EUV adjustment
#     rho_euv = rho_thermo * (F10.7 / 100)**0.2
    
#     # Dust storm scaling
#     rho_dust = rho_euv * (1 + 0.15 * dust_tau)
    
#     return rho_dust

def mars_beta_angle(season: str, inclination_deg: float = 93.0, orbit_altitude_km: float = 200.0) -> float:
    """
    Compute the beta angle for a satellite orbiting Mars, based on season, inclination, and altitude.

    Parameters:
        season (str): 'spring', 'summer', 'autumn', or 'winter'
        inclination_deg (float): Orbital inclination in degrees
        orbit_altitude_km (float): Orbital altitude above Mars in km (not used in beta calc, but included for completeness)

    Returns:
        beta_deg (float): Beta angle in degrees
    """
    # Mars constants
    epsilon_deg = 25.19  # Mars axial tilt in degrees

    # Map seasons to Mars solar longitude (Ls)
    season_to_Ls = {
        "spring": 0,
        "summer": 90,
        "autumn": 180,
        "fall": 180,  # alias
        "winter": 270
    }

    season = season.lower()
    if season not in season_to_Ls:
        raise ValueError("Season must be one of: spring, summer, autumn (or fall), winter.")

    # Compute solar declination δ from Ls
    Ls = season_to_Ls[season]
    delta_deg = epsilon_deg * np.sin(np.radians(Ls))

    # Compute beta angle
    sin_beta = np.sin(np.radians(delta_deg)) / np.cos(np.radians(inclination_deg))
    sin_beta = np.clip(sin_beta, -1, 1)  # Prevent math domain errors
    beta_rad = np.arcsin(sin_beta)
    beta_deg = np.degrees(beta_rad)

    return beta_deg

if __name__ == "__main__":
    # Example usage
    h = 200  # altitude in km
    eclipse_time, eclipse_percent = calc_2d_elipse(h)
    print(f"Eclipse time: {eclipse_time:.2f} s")
    print(f"Eclipse percent: {eclipse_percent:.2f} %")

    # Example beta angle calculation
    beta = mars_beta_angle(season="autumn", inclination_deg=93, orbit_altitude_km=200)
    print(f"Beta angle: {beta}°")