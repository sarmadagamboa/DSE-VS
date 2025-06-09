import numpy as np
import matplotlib.pyplot as plt

def calc_drag_strategy(start_h, margin, mission_lifetime_years=1.9):
    mu_mars = 42828.376645  # km^3 / s^2
    R_M = 3396.2  # km
    sol_minutes = 24.623 * 60  # min
    a = R_M + start_h  # km

    h_max = start_h + margin
    h_min = start_h - margin

    delta_h_per_orbit = 8.37 / 1000  # km

    a_max = R_M + h_max
    a_min = R_M + h_min

    period = 2 * np.pi * np.sqrt(a**3 / mu_mars)
    period_minutes = period / 60

    delta_h_total = h_max - h_min
    num_orbits = delta_h_total / delta_h_per_orbit

    total_minutes = num_orbits * period_minutes
    total_sols = total_minutes / sol_minutes

    v_circ_min = np.sqrt(mu_mars / a_min)
    v_circ_max = np.sqrt(mu_mars / a_max)
    
    v_transfer_periapsis = np.sqrt(mu_mars * (2 / a_min - 1 / ((a_min + a_max) / 2)))
    v_transfer_apoapsis  = np.sqrt(mu_mars * (2 / a_max - 1 / ((a_min + a_max) / 2)))

    delta_v1 = v_transfer_periapsis - v_circ_min
    delta_v2 = v_circ_max - v_transfer_apoapsis

    delta_v_total = (delta_v1 + delta_v2) * 1000

    mission_lifetime_sols = mission_lifetime_years * 365.25 / 1.0275
    num_burns = mission_lifetime_sols / total_sols
    lifetime_delta_v = delta_v_total * num_burns

    print(f"Margin: {margin:.3f} km | ΔV: {delta_v_total:.3f} m/s | "
          f"Δt between burns: {total_sols:.2f} sols | "
          f"Num burns: {num_burns:.1f} | Lifetime ΔV: {lifetime_delta_v:.2f} m/s")

    return delta_v_total, total_sols, lifetime_delta_v

start_h = 212.48 

margins = np.linspace(0.01, 1.0, 50) 

delta_v_list = []
time_sols_list = []

for margin in margins:
    delta_v, total_sols, lifetime_delta_v = calc_drag_strategy(start_h, margin)
    delta_v_list.append(delta_v)
    time_sols_list.append(total_sols)
    
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(margins*1000, delta_v_list, marker='o')
plt.xlabel("Margin (m)")
plt.ylabel("Delta-V (m/s)")
plt.title("Delta-V vs Margin")
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(margins*1000, time_sols_list, marker='s', color='orange')
plt.xlabel("Margin (m)")
plt.ylabel("Time between burns (sols)")
plt.title("Time vs Margin")
plt.grid()

plt.tight_layout()
plt.show()

def precession_rate(a, inc_deg):

    mu_mars = 42828.376645  # km^3 / s^2
    R_M = 3396.2  # km
    J2 = 0.001958705252674          # constant of oblateness, J2

    i = np.radians(inc_deg)

    return -(3/2) * J2 * (R_M**2 / a**2) * np.sqrt(mu_mars / a**3) * np.cos(i)

def precession_diff_percent(start_h, margin):

    # Constants
    mu_mars = 42828.376645  # km^3/s^2
    R_M = 3396.2  # km
    J2 = 0.001958705252674          # constant of oblateness, J2
    inc_deg = 92.4  # degrees

    a_nominal = R_M + start_h
    a_plus = R_M + start_h + margin
    a_minus = R_M + start_h - margin

    rate_nominal = precession_rate(a_nominal, inc_deg)
    rate_plus = precession_rate(a_plus, inc_deg)
    rate_minus = precession_rate(a_minus, inc_deg)

    diff_plus = (rate_plus - rate_nominal) / rate_nominal * 100
    diff_minus = (rate_minus - rate_nominal) / rate_nominal * 100

    return diff_plus, diff_minus

def accumulated_raan_drift(delta_percent, mission_years=1.9):

    delta_frac = delta_percent / 100.0
    
    sun_sync_rate_deg_day = 360 / 687.0 
    
    mission_days = mission_years * 365.25
    
    delta_RAAN_deg = delta_frac * sun_sync_rate_deg_day * mission_days
    
    return delta_RAAN_deg

for margin in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    diff_plus, diff_minus = precession_diff_percent(212.48, margin)
    drift_plus = accumulated_raan_drift(diff_plus)
    drift_minus = accumulated_raan_drift(diff_minus)
    print(f"Margin ±{margin*1000:.0f} m: {drift_plus:.4f} deg, {drift_minus:.4f} deg")