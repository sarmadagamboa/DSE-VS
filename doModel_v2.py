import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def plot_sst():
    # Parameters
    h = np.array([100000, 150000, 200000, 250000, 300000, 350000])  # satellite altitudes (m)
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    d_fixed = 200000  # satellite separation (m)
    f = 1  # observation frequency (Hz) 1 measurement per second

    t_5_days = 5 * 24 * 3600
    t_30_days = 30 * 24 * 3600

    epsilon = 1e-8
    colors = ['b', 'g', 'r', 'purple', 'orange', 'cyan']
    l_range = np.arange(1, 401)

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    durations = [(t_5_days, "5 Days", axs[0]), (t_30_days, "30 Days", axs[1])]

    for t, label, ax in durations:
        N = f * t
        for i, height in enumerate(h):
            current_r = R + height
            current_n = n[i]

            gamma = 2* np.arctan((d_fixed/2)/current_r)

            #Fl = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2)) * (1 / np.sqrt(l_range)) #smoothened 

            #Fl = l_range / np.sqrt(2 + 4 * l_range + 4 * l_range**2) * (1 / np.abs(np.sin(l_range * gamma / 2)))
            #sigma = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * (current_r / R)**l_range * Fl

            #try stopping curve at first anomaly and then integrate
            # Compute raw F_l with possible anomalies

            Fl_raw = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
            with np.errstate(divide='ignore', invalid='ignore'):
                Fl_raw *= (1 / np.abs(np.sin(l_range * gamma / 2)))

            # --- Robust Spike Detection using Gradient and Peak Detection ---
            log_Fl_raw = np.log10(Fl_raw + 1e-20)  # avoid log(0)
            gradient = np.gradient(log_Fl_raw)
            peaks, _ = find_peaks(gradient, height=0.5)  # threshold may be tuned

            # Fallback if no spike found
            if len(peaks) > 0:
                spike_idx = peaks[0]
            else:
                spike_idx = len(Fl_raw)

            # --- Fit Tail After Spike ---
            Fl_valid = Fl_raw[:spike_idx + 1]
            l_valid = l_range[:spike_idx + 1]
            l_tail = l_range[spike_idx + 1:]

            if len(l_tail) > 0:
                fit_l = l_valid[-5:]
                fit_fl = Fl_valid[-5:]
                coeffs = np.polyfit(fit_l, np.log(fit_fl), 1)
                Fl_tail = np.exp(coeffs[1] + coeffs[0] * l_tail)
                Fl_final = np.concatenate([Fl_valid, Fl_tail])
            else:
                Fl_final = Fl_valid

            # --- Compute Sigma ---
            sigma_raw = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * (current_r / R)**l_range * Fl_raw
            sigma = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * (current_r / R)**l_range * Fl_final

            # --- Plotting ---
            ax.plot(l_range, sigma, color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')
            ax.plot(l_range, sigma_raw, color=colors[i], linestyle=':', alpha=0.6)

            #ax.plot(l_range, sigma, color=colors[i], linewidth=2, label=f'{height/1000:.0f} km')

        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.tight_layout()
    plt.show()

# Run the function
plot_sst()
