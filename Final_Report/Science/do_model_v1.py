import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import matplotlib.cm as cm

def plot_sst_qt():
    # Parameters
    #h = np.array([100000, 150000, 200000, 250000, 300000])  # satellite altitudes (m)
    h = np.array([153800,177000, 188000, 212000])

    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    d_fixed = 60000  # satellite separation (m) #220km GRACE-FO
    f = 1  # observation frequency (Hz)

    SOL_DURATION = 88775.244  # seconds in one Martian sol

    t_5_sols = 5 * SOL_DURATION
    t_30_sols = 30 * SOL_DURATION
    #t_1_sols = 668.6 * SOL_DURATION
    t_3_years = 3 * 365 * 24 * 60 *60

    epsilon = 1e-10
    colors = ['b', 'g', 'r', 'purple', 'orange', 'cyan']
    l_range = np.arange(1, 401)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    durations = [(t_5_sols, "5 Sols", axs[0],50), (t_30_sols, "30 Sols", axs[1],90), (t_3_years, "3 Earth Years ", axs[2],90)]

    for t, label, ax, cutoff in durations:
        N = f * t
        for i, height in enumerate(h):
            current_r = R + height
            current_n = n[i]

            gamma = 2 * np.arctan((d_fixed / 2) / current_r)

            # Just compute the altitude lines directly
            with np.errstate(divide='ignore', invalid='ignore'):
                Fl = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
                Fl *= (1 / np.abs(np.sin(l_range * gamma / 2) + 1e-10))
            
            # Calculate sigma
            sigma = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * ((current_r / R)**l_range) * Fl

            ax.plot(l_range, sigma, color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')

        # Plot the reference signal power law
        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-18, top=1e-4)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity SST measurement with LRI+CAI')
    plt.tight_layout()
    plt.show()

plot_sst_qt()