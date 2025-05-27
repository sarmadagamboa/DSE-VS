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

#