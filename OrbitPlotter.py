from astropy import units as u
from astropy.constants import G, M_earth, R_earth
import numpy as np
from math import floor, pi
import matplotlib.pyplot as plt


# Add these lines manually to define Mars constants
M_mars = 6.4169e23   # Mars mass in kg
R_mars = 3396.2e3        # Mars radius (mean)
# Constants
mu_mars = 6.67e-11 * M_mars  # km^3/s^2 (Mars)
Tm = 88642 #Mars sidereal day length in seconds
mars_J2 = 1.9555e-3



def a_to_T(a, mu):
    """
    Calculate the period of an orbit given its semi-major axis.

    Parameters:
    a: Semi-major axis in m.

    Returns:
    float: Orbital period in seconds.
    """
    return np.sqrt((4 * a**3 * np.pi**2) / mu)


def repeat_orbit_e(a, e, k):
    """
    Calculate possible repeat orbits for a given semi-major axis and eccentricity.

    Parameters:
    a (float): Semi-major axis in m
    e (float): Eccentricity of the orbit
    k (int): Number of days between repeats

    Returns:
    possible inclinations (list): List of possible inclinations in degrees
    """
    
    T = a_to_T(a, mu_mars)

    rep_time = k * Tm
    
    j1 = np.floor(rep_time / T)
    print(j1)
    j2 = j1 + 1
    print(j2)

    i1p = np.arccos((+(k*2*np.pi/j1)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i1p)
    i1m = np.arccos((-(k*2*np.pi/j1)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i1m)
    i2p = np.arccos((+(k*2*np.pi/j2)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i2p)
    i2m = np.arccos((-(k*2*np.pi/j2)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i2m)

    plt.plot(e, np.degrees(i1p), label='i1p')
    plt.plot(e, np.degrees(i1m), label='i1m')
    plt.plot(e, np.degrees(i2p), label='i2p')
    plt.plot(e, np.degrees(i2m), label='i2m')
    plt.title('Inclination vs Eccentricity for Repeat Orbits')
    plt.xlabel('Eccentricity')
    plt.ylabel('Inclination (degrees)')
    plt.legend()
    plt.show()


def repeat_orbit_a(a, e, k):
    """
    Calculate possible repeat orbits for a given semi-major axis and eccentricity.

    Parameters:
    a (float): Semi-major axis in m
    e (float): Eccentricity of the orbit
    k (int): Number of days between repeats

    Returns:
    possible inclinations (list): List of possible inclinations in degrees
    """
    
    T = a_to_T(a, mu_mars)

    rep_time = k * Tm
    
    j1 = np.floor(rep_time / T)
    print(j1)
    j2 = j1 + 1
    print(j2)

    i1p = np.arccos(np.multiply((+(k*2*np.pi/j1)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm), ((np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))))
    print(i1p)
    i1m = np.arccos((-(k*2*np.pi/j1)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i1m)
    i2p = np.arccos((+(k*2*np.pi/j2)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i2p)
    i2m = np.arccos((-(k*2*np.pi/j2)-2*np.pi*(2*np.pi*np.sqrt(mu_mars/np.power(a, 3)))/Tm)*(np.power(a, 2)*(1-np.power(e, 2))**2)/(3*np.pi*mars_J2*R_mars**2))
    print(i2m)

    plt.plot(a, np.degrees(i1p), label='i1p')
    plt.plot(a, np.degrees(i1m), label='i1m')
    plt.plot(a, np.degrees(i2p), label='i2p')
    plt.plot(a, np.degrees(i2m), label='i2m')
    plt.title('Inclination vs Eccentricity for Repeat Orbits')
    plt.xlabel('Eccentricity')
    plt.ylabel('Inclination (degrees)')
    plt.legend()
    plt.show()


#repeat_orbit(R_mars+220e3, 0.2, 1)


e = np.linspace(0, 0.2, 50)
a = 300e3
k = 3
#repeat_orbit_e(a, e, k)

e = 0
a = np.linspace(R_mars+220e3, R_mars+220e3+200e3, 50)
k = 3
repeat_orbit_a(a, e, k)