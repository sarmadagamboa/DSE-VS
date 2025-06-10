import numpy as np
import matplotlib.pyplot as plt

# === Parameters ===
battery_capacity_Wh = 1477   # Wh
initial_soc = 1            

orbital_period_min = 109
sunlit_fraction = 0.976
orbital_period_s = orbital_period_min * 60

P_draw_sun = 820.82   # W
P_draw_ecl = 371.02   # W
ecl_fraction = 1 - sunlit_fraction
P_gen = 980

# Peak‐draw event
peak_period_s       = 57 * 3600
peak_extra_W        = 1025
peak_event_duration = 20 * 60

# Simulation grid
total_sim_hours = 60
dt = 60  # seconds
t = np.arange(0, total_sim_hours * 3600, dt)

# Compute net power profile
net_power = np.zeros_like(t, dtype=float)
for i, tt in enumerate(t):
    sec_in_orbit = tt % orbital_period_s
    # Charge during sunlight, discharge during eclipse
    if sec_in_orbit < sunlit_fraction * orbital_period_s:
        net_power[i] = P_gen - P_draw_sun
    else:
        net_power[i] = -P_draw_ecl
    # Add peak draw
    if (tt % peak_period_s) < peak_event_duration:
        net_power[i] -= peak_extra_W

energy = np.zeros_like(t, dtype=float)
soc = np.zeros_like(t, dtype=float)

energy[0] = battery_capacity_Wh * initial_soc
soc[0]    = initial_soc

for i in range(1, len(t)):
    energy[i] = energy[i-1] + net_power[i] * dt / 3600.0
    energy[i] = np.clip(energy[i], 0, battery_capacity_Wh)
    soc[i]    = energy[i] / battery_capacity_Wh

# Plot
time_h = t / 3600.0
plt.figure(figsize=(10, 4))
plt.plot(time_h, soc * 100, label='State of Charge [%]')
plt.xlabel('Time [hours]')
plt.ylabel('SoC [%]')
plt.title(f'Battery State of Charge Over {total_sim_hours:.1f} h')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
