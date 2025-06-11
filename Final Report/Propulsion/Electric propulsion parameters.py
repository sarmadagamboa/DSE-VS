import numpy as np
import matplotlib.pyplot as plt

def compute_average_thrust(m_initial, deltav, duration_minutes):
    duration_s = duration_minutes * 60
    return m_initial * deltav / duration_s

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


# Parameters
m_initial = 650.0  # kg
start_h = 212.48  # km

# Define range of Delta-Vs and Burn Durations
delta_vs = np.linspace(0.01, 0.6, 200)  # m/s
durations = np.linspace(1, 180, 300)     # minutes
DELTA_VS, DURATIONS_1 = np.meshgrid(delta_vs, durations)
# Compute thrusts
THRUST_1 = compute_average_thrust(m_initial, DELTA_VS, DURATIONS_1)
# Mask values above 0.4 N (400 mN)
THRUST_1_MASKED = np.ma.masked_where(THRUST_1 > 0.250, THRUST_1)

#  Define range of Margin and Burn Durations
margins = np.linspace(0.01, 0.8, 200)  # km
MARGINS, DURATIONS_2 = np.meshgrid(margins, durations)
DELTA_V_FROM_MARGIN = np.vectorize(lambda m: calc_drag_strategy(start_h, m)[0])(MARGINS)
THRUST_2 = compute_average_thrust(m_initial, DELTA_V_FROM_MARGIN, DURATIONS_2)
THRUST_2_MASKED = np.ma.masked_where(THRUST_2 > 0.250, THRUST_2)

fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Delta-V vs Duration
contour1 = axs[0].contourf(DELTA_VS, DURATIONS_1, THRUST_1_MASKED * 1000, levels=50, cmap='plasma')
cbar1 = fig.colorbar(contour1, ax=axs[0])
cbar1.set_label("Thrust (mN)")
axs[0].set_title("Thrust vs Delta-V and Burn Duration")
axs[0].set_xlabel("Delta-V (m/s)")
axs[0].set_ylabel("Burn Duration (min)")
axs[0].grid(True)

# Plot 2: Margin vs Duration
contour2 = axs[1].contourf(MARGINS * 1000, DURATIONS_2, THRUST_2_MASKED * 1000, levels=50, cmap='plasma')
cbar2 = fig.colorbar(contour2, ax=axs[1])
cbar2.set_label("Thrust (mN)")
axs[1].set_title("Thrust vs Margin and Burn Duration")
axs[1].set_xlabel("Margin (m)")
axs[1].set_ylabel("Burn Duration (min)")
axs[1].grid(True)

plt.tight_layout()
plt.show()
