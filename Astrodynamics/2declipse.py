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