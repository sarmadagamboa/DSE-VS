import numpy as np
import matplotlib.pyplot as plt

def calc_drag_strategy(start_h, margin):
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

    v_circ_max = np.sqrt(mu_mars / a_max)
    v_circ_min = np.sqrt(mu_mars / a_min)
    delta_v = abs(v_circ_max - v_circ_min) * 1000  # m/s

    return delta_v, total_sols

start_h = 212.48 

margins = np.linspace(0.01, 1.0, 50) 

delta_v_list = []
time_sols_list = []

for margin in margins:
    delta_v, total_sols = calc_drag_strategy(start_h, margin)
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