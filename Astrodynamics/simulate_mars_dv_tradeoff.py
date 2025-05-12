
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import G, M_earth, R_earth
from astropy import units as u

# === Reuse constants and mission functions from the final user code ===
# (Assumes the functions are already defined in this context)

# Import your main and function definitions here if in a separate file/module
# from your_module import main

# Constants (reuse from user definitions)
M_mars = 6.4171e23 * u.kg
R_mars = 3389.5 * u.km
d_MS = 228e6
d_ES = 149e6
mu_earth = (G * M_earth).to(u.km**3 / u.s**2).value
mu_mars = (G * M_mars).to(u.km**3 / u.s**2).value
r_earth = R_earth.to(u.km).value
r_mars = R_mars.to(u.km).value
mu_sun = 1.32712440018e11

# Import compute and main functions from final user-provided script
from mars_delta_v_budget_functions import main  # Update this path if needed

# Define altitudes
leo_altitudes = np.arange(150, 501, 25)
mars_altitudes = np.arange(200, 501, 25)

# Define condition set
conditions = [
    {"use_aerobraking": False, "inclination_midcourse": False},
    {"use_aerobraking": True,  "inclination_midcourse": False},
    {"use_aerobraking": False, "inclination_midcourse": True},
    {"use_aerobraking": True,  "inclination_midcourse": True},
]

results = {}

# Simulation loop
for cond in conditions:
    Z = np.zeros((len(leo_altitudes), len(mars_altitudes)))
    for i, leo_alt in enumerate(leo_altitudes):
        for j, mars_alt in enumerate(mars_altitudes):
            try:
                import io
                import contextlib
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    main(
                        use_aerobraking=cond["use_aerobraking"],
                        inclination_midcourse=cond["inclination_midcourse"],
                        leo_alt=leo_alt,
                        mars_orbit_alt=mars_alt
                    )
                output = buf.getvalue()
                for line in output.splitlines():
                    if "Total Mission ΔV" in line:
                        dv = float(line.split('=')[1].split('km/s')[0].strip())
                        Z[i, j] = dv
                        break
            except Exception as e:
                Z[i, j] = np.nan
    results[(cond["use_aerobraking"], cond["inclination_midcourse"])] = Z

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
titles = [
    "No Aerobraking, No Midcourse",
    "Aerobraking, No Midcourse",
    "No Aerobraking, Midcourse",
    "Aerobraking, Midcourse"
]

for ax, ((ab, mc), Z), title in zip(axs.flat, results.items(), titles):
    im = ax.imshow(
        Z.T,
        extent=[leo_altitudes.min(), leo_altitudes.max(), mars_altitudes.min(), mars_altitudes.max()],
        origin='lower',
        aspect='auto',
        cmap='viridis'
    )
    ax.set_title(title)
    ax.set_xlabel("LEO Altitude (km)")
    ax.set_ylabel("Mars Orbit Altitude (km)")
    fig.colorbar(im, ax=ax, label='Total ΔV (km/s)')

plt.suptitle("ΔV Budget vs LEO and Mars Orbit Altitudes (All Scenarios)", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
