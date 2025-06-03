import numpy as np

def compute_orbital_period(a_km, mu):
    """
    Compute orbital period in days from semi-major axis and gravitational parameter.
    a_km: semi-major axis in km
    mu: gravitational parameter of central body in km^3/s^2
    """
    period_seconds = 2 * np.pi * np.sqrt(a_km**3 / mu)
    period_days = period_seconds / (3600 * 24)
    return period_days

def compute_synodic_period(T_inner, T_outer):
    return 1 / abs(1 / T_inner - 1 / T_outer)

mu_sun = 1.32712440018e11  # km^3/s^2, standard gravitational parameter of the Sun

a_earth = 149.6e6  # km
a_mars = 227.9e6   # km

T_earth = compute_orbital_period(a_earth, mu_sun)
T_mars = compute_orbital_period(a_mars, mu_sun)

S = compute_synodic_period(T_earth, T_mars)

print(f"Earth orbital period: {T_earth:.2f} days")
print(f"Mars orbital period: {T_mars:.2f} days")
print(f"Synodic period Earth-Mars: {S:.2f} days (~{S/365.25:.2f} years)")
