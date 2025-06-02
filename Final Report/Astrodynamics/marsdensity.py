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

    mask = Ls == 270
    dens_270 = density[mask]
    height_270 = height[mask]

    H = np.arange(80, 245, 5)  # 80:5:240
    dmax = np.zeros_like(H, dtype=float)
    dmin = np.zeros_like(H, dtype=float)

    for i, h in enumerate(H):
        dens_at_h = dens_270[height_270 == h]
        if len(dens_at_h) > 0:
            dmax[i] = np.max(dens_at_h)
            dmin[i] = np.min(dens_at_h)
        else:
            dmax[i] = np.nan
            dmin[i] = np.nan

    plt.figure(figsize=(8, 6))
    plt.plot(np.log10(dmax), H, 'k', linewidth=4, label='max density')
    plt.scatter(np.log10(dens_270), height_270, color='r', marker='.', label='range of density')
    plt.plot(np.log10(dmin), H, 'b', linewidth=4, label='min density')
    plt.xlabel('density (log10 kg/m³)')
    plt.ylabel('height (km)')
    if year == 2:
        plt.title("Worst-Case Mars Atmospheric Density vs Altitude\n(Ls = 270°, High Solar, TES Year 2)")
    elif year == 1:
        plt.title("Worst-Case Mars Atmospheric Density vs Altitude\n(Ls = 270°, High Solar, TES Year 1)")
    plt.legend()
    plt.grid(True)
    plt.show()

plot_mars_density(year=2)