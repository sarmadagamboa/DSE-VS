import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def plot_gravity_aquifers_single_plot():
    # Constants
    G = 6.67430e-11  # gravitational constant, m^3/kg/s^2
    rho = 997 - 3050  # density contrast (water - crust), kg/m^3
    pi = np.pi

    # Mars radius
    R = 3389.5e3  # m

    # Depths at different Mars regions
    Depths = {
        'Equator (4 km depth)': 4000,  
        'Mid Latitude (8 km depth)': 8000,            
        'Polar (18 km depth)': 18000,                              
    }

    # Aquifer thickness from 10 m to 3000 m
    h_range = np.linspace(10, 3000, 500)

    # Define colors manually
    colors = ['firebrick', 'darkorange', 'goldenrod']

    # Create plot
    plt.figure(figsize=(10, 6))

    for i, (label, D) in enumerate(Depths.items()):
        r = R + D
        delta_g = 2 * pi * rho * h_range * G * (R / r)  # in m/s²
        delta_g_ugal = delta_g * 1e6  # convert to μGal

        plt.plot(h_range, delta_g_ugal, label=label, color=colors[i], linewidth=2)

    # Labels and title
    plt.xlabel('Aquifer Thickness (m)', fontsize=12)
    plt.ylabel('Gravity Anomaly (μGal)', fontsize=12)
    plt.title('Mars Gravity Anomaly vs Aquifer Thickness\nby Regional Aquifer Depth', 
              fontsize=14, fontweight='bold')

    # Legend and grid
    plt.legend(title='Region / Depth', fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.xlim(0, 3100)
    

    # Final layout
    plt.tight_layout()

    # Console output for reference values
    print("Gravity Anomaly Analysis for Mars Aquifers")
    print("=" * 50)
    print(f"Density contrast (ρ): {rho} kg/m³")
    print(f"Mars radius (R): {R/1000:.1f} km\n")

    for label, D in Depths.items():
        r = R + D
        delta_g_ref = 2 * pi * rho * 1000 * G * (R / r) * 1e6  # μGal
        print(f"{label}:")
        print(f"  Gravity anomaly for 1000m aquifer: {delta_g_ref:.2f} μGal\n")

    
    # Show plot
    plt.show()

def plot_gravity_anomaly_difference():
    # Constants
    G = 6.67430e-11  # gravitational constant
    rho = 997 - 3050  # density contrast (water - crust)
    pi = np.pi
    R = 3389.5e3  # Mars radius (m)

    # Aquifer thickness range
    h_range = np.linspace(10, 3000, 500)

    # Region depths
    Depths = {
        'Equator': 4000,
        'Mid Latitude': 8000,
        'Polar': 18000,
    }

    # Choose Equator as baseline
    ref_D = Depths['Equator']
    ref_r = R + ref_D
    ref_delta_g = 2 * pi * rho * h_range * G * (R / ref_r) * 1e6  # in μGal

    # Plot setup
    plt.figure(figsize=(10, 6))
    colors = {'Mid Latitude': 'darkorange', 'Polar': 'firebrick'}

    for region, D in Depths.items():
        if region == 'Equator':
            continue  # skip baseline

        r = R + D
        delta_g = 2 * pi * rho * h_range * G * (R / r) * 1e6  # in μGal
        delta_diff = delta_g - ref_delta_g  # Δ from Equator

        plt.plot(h_range, delta_diff, label=f'{region} - Equator', 
                 color=colors[region], linewidth=2)

    # Plot formatting
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.xlabel('Aquifer Thickness (m)', fontsize=12)
    plt.ylabel('Δ Gravity Anomaly (μGal)', fontsize=12)
    plt.title('Difference in Mars Gravity Anomaly\nCompared to Equator Region', 
              fontsize=14, fontweight='bold')
    plt.legend(title='Δ from Equator')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.xlim(0, 3100)

    # Show plot
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def plot_equator_and_deltas():
    # ---------- constants ----------
    G   = 6.67430e-11                    # m³ kg⁻¹ s⁻²
    rho = 997 - 3050                     # kg m⁻³  (water – crust)
    pi  = np.pi
    R   = 3389.5e3                       # m  (Mars radius)

    # ---------- aquifer thickness range ----------
    h_range = np.linspace(10, 3000, 500)        # 10 m … 3 km

    # ---------- region depths ----------
    depths = {
        'Equator':      4000,   # m
        'Mid Latitude': 8000,
        'Polar':        18000,
    }

    # ---------- baseline (Equator) ----------
    r_eq  = R + depths['Equator']
    g_eq  = 2 * pi * rho * h_range * G * (R / r_eq) * 1e6     # μGal

    # ---------- make two‐panel figure ----------
    fig, (ax_abs, ax_delta) = plt.subplots(
        1, 2, figsize=(13, 5), sharex=True, gridspec_kw={'wspace': 0.28}
    )

    # === LEFT: absolute anomaly (Equator only) ===
    ax_abs.plot(h_range, g_eq, color='firebrick', lw=2)
    ax_abs.set_title('Equator Aquifer\nAbsolute Gravity Anomaly', fontweight='bold')
    ax_abs.set_xlabel('Aquifer Thickness (m)')
    ax_abs.set_ylabel('Gravity Anomaly (μGal)')
    ax_abs.grid(ls='--', alpha=0.5)
    ax_abs.set_xlim(0, 3100)

    # === RIGHT: Δ relative to Equator ===
    color_map = {'Mid Latitude': 'darkorange', 'Polar': 'goldenrod'}
    for region, D in depths.items():
        if region == 'Equator':
            continue                      # skip baseline
        r     = R + D
        g_reg = 2 * pi * rho * h_range * G * (R / r) * 1e6
        delta = g_reg - g_eq              # μGal difference
        ax_delta.plot(h_range, delta,
                      label=f'{region} – Equator',
                      color=color_map[region], lw=2)

    ax_delta.axhline(0, color='grey', ls='--', lw=1)
    ax_delta.set_title('Δ Gravity Anomaly\nvs Equator Baseline', fontweight='bold')
    ax_delta.set_xlabel('Aquifer Thickness (m)')
    ax_delta.set_ylabel('Δ Gravity Anomaly (μGal)')
    ax_delta.grid(ls='--', alpha=0.5)
    ax_delta.set_xlim(0, 3100)
    ax_delta.legend(frameon=False, fontsize=9)

    
    fig.suptitle('Mars Aquifer Gravity-Anomaly Comparison', fontsize=15, fontweight='bold', y=1.05)
    fig.tight_layout()
    plt.show()




# Run 
plot_gravity_aquifers_single_plot()
#plot_gravity_anomaly_difference()
plot_equator_and_deltas()
