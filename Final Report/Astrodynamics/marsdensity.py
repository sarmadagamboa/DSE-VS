import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_mars_density(year=2):
    if year == 2:
        data = np.loadtxt(r'Final Report\Astrodynamics\tpdhsy21.txt', skiprows=1) 
    elif year == 1:
        data = np.loadtxt(r'Final Report\Astrodynamics\tpdhsy11.txt', skiprows=1) 

    Ls = data[:, 0]
    height = data[:, 1]
    density = data[:, 13]

    H = np.arange(80, 300, 5)

    dmax_all = np.zeros_like(H, dtype=float)
    dmin_all = np.zeros_like(H, dtype=float)
    for i, h in enumerate(H):
        dens_at_h = density[height == h]
        if len(dens_at_h) > 0:
            dmax_all[i] = np.max(dens_at_h)
            dmin_all[i] = np.min(dens_at_h)
        else:
            dmax_all[i] = np.nan
            dmin_all[i] = np.nan

    mask_270 = Ls == 270
    dens_270 = density[mask_270]
    height_270 = height[mask_270]

    dmax_270 = np.zeros_like(H, dtype=float)
    dmin_270 = np.zeros_like(H, dtype=float)
    for i, h in enumerate(H):
        dens_at_h_270 = dens_270[height_270 == h]
        if len(dens_at_h_270) > 0:
            dmax_270[i] = np.max(dens_at_h_270)
            dmin_270[i] = np.min(dens_at_h_270)
        else:
            dmax_270[i] = np.nan
            dmin_270[i] = np.nan

    plt.figure(figsize=(12, 6))

    plt.scatter(np.log10(density), height, color='r', s=5, label='Density range (all Ls)')

    plt.plot(np.log10(dmax_all), H, 'k-', linewidth=3, label='Max density (all Ls)')
    plt.plot(np.log10(dmin_all), H, 'b-', linewidth=3, label='Min density (all Ls)')

    plt.plot(np.log10(dmax_270), H, 'k--', linewidth=3, label='Max density (Ls=270°)')
    plt.plot(np.log10(dmin_270), H, 'b--', linewidth=3, label='Min density (Ls=270°)')

    plt.xlabel('density (log10 kg/m³)')
    plt.ylabel('height (km)')
    if year == 2:
        plt.title("Mars Atmospheric Density Envelope (TES Year 2)")
    elif year == 1:
        plt.title("Mars Atmospheric Density Envelope (TES Year 1)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Return max densities for further use
    return {
        'height': H,
        'max_density': dmax_all
    }

# Example usage
results = plot_mars_density(year=2)
print("Heights (km):", results['height'])
print("Max Densities (kg/m^3):", results['max_density'])
