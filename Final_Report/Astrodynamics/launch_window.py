import numpy as np
from datetime import datetime, timedelta

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
    """
    Compute synodic period in days between two orbital periods (T_inner < T_outer).
    """
    return 1 / abs(1 / T_inner - 1 / T_outer)

def compute_phase_angle(a_inner, a_outer):
    """
    Compute ideal phase angle (degrees) for Hohmann transfer between orbits with semi-major axes a_inner and a_outer
    """
    return 180 * (1 - np.sqrt(((a_inner+a_outer)**3)/(8*(a_outer)**3)))

def true_phase_angle(t, T_inner, T_outer):
    """
    Compute Earth-Mars phase angle at time t (days) since reference zero point.
    Returns angle in degrees between 0 and 360.
    """
    mean_motion_inner = 360 / T_inner  
    mean_motion_outer = 360 / T_outer

    # Phase angle is difference in angular positions
    angle = (mean_motion_outer - mean_motion_inner) * t + 173.4
    return angle % 360  # Normalize to [0, 360) degrees

def compute_launch_windows(start_date, synodic_period, T_earth, T_mars, ideal_phase_angle, tolerance_deg=1, step_days=1):
    """
    Compute launch windows over one synodic period starting from start_date (datetime).
    
    Returns a list of datetime objects for valid launch windows.
    """
    windows = []
    for day_offset in range(int(synodic_period) + 1):
        t = day_offset  # days since start_date
        phase = true_phase_angle(t, T_earth, T_mars)
        # Check if phase angle is within tolerance of ideal phase angle or its complementary angle
        if (abs(phase - ideal_phase_angle) <= tolerance_deg) or (abs(phase - (360 - ideal_phase_angle)) <= tolerance_deg):
            launch_day = start_date + timedelta(days=day_offset)
            windows.append((launch_day, phase))
    return windows

# Constants
mu_sun = 1.32712440018e11  # km^3/s^2

a_earth = 149.6e6  # km
a_mars = 227.9e6   # km

# Calculate orbital periods
T_earth = compute_orbital_period(a_earth, mu_sun)
T_mars = compute_orbital_period(a_mars, mu_sun)

# Synodic period in days
S = compute_synodic_period(T_earth, T_mars)

# Ideal phase angle for transfer (degrees)
ideal_phase = compute_phase_angle(a_earth, a_mars)
print(f"Earth orbital period: {T_earth:.2f} days")
print(f"Mars orbital period: {T_mars:.2f} days")
print(f"Synodic period Earth-Mars: {S:.2f} days (~{S/365.25:.2f} years)")
print(f"Ideal phase angle for transfer: {ideal_phase:.2f} degrees")

# Compute launch windows starting today
start = datetime(2041, 1, 1)
launch_windows = compute_launch_windows(start, S, T_earth, T_mars, ideal_phase, tolerance_deg=1, step_days=1)

print("\nLaunch windows (within ±1 degrees of ideal phase angle):")
for win, phase in launch_windows:
    print(f"Date: {win.strftime('%Y-%m-%d')}, Phase angle: {phase:.2f}°")
    print(f'Arrival Date: {(win + timedelta(days=258.88)).strftime("%Y-%m-%d")}, Phase angle: {true_phase_angle(S, T_earth, T_mars):.2f}°')
