import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import matplotlib.cm as cm

def plot_sst():
    # Parameters
    h = np.array([100000, 150000, 200000, 250000, 300000])  # satellite altitudes (m)
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    d_fixed = 220000  # satellite separation (m) #220km GRACE-FO
    f = 1  # observation frequency (Hz)

    SOL_DURATION = 88775.244  # seconds in one Martian sol

    t_5_sols = 5 * SOL_DURATION
    t_30_sols = 30 * SOL_DURATION
    t_1_sols = 668.6 * SOL_DURATION

    epsilon = 1e-8
    colors = ['b', 'g', 'r', 'purple', 'orange', 'cyan']
    l_range = np.arange(1, 401)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    durations = [(t_5_sols, "5 Sols", axs[0],50), (t_30_sols, "30 Sols", axs[1],90), (t_1_sols, "1 Martian Year ", axs[2],90)]

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
        ax.set_ylim(bottom=1e-16, top=1e-4)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity SST measurement with LRI+ACC')
    plt.tight_layout()
    plt.show()



def plot_sst_v2():
    # Parameters
    h = np.array([100000, 150000, 200000, 250000, 300000, 350000])  # satellite altitudes (m)
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    d_fixed = 200000  # satellite separation (m)
    f = 1  # observation frequency (Hz)

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

            gamma = 2 * np.arctan((d_fixed / 2) / current_r)

            # Compute raw F_l with possible anomalies
            with np.errstate(divide='ignore', invalid='ignore'):
                Fl_raw = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
                Fl_raw *= (1 / np.abs(np.sin(l_range * gamma / 2) + 1e-10))

                def extrapolate_tail(x, y, start_point):
                    fit_x = x[:start_point][-5:]
                    fit_y = y[:start_point][-5:]
                    log_y = np.log(fit_y + 1e-10)
                    coeffs = np.polyfit(fit_x, log_y, 1)
                    tail_x = x[start_point:]
                    tail_y = np.exp(coeffs[1] + coeffs[0] * tail_x)
                    return tail_y

            log_Fl_raw = np.log10(Fl_raw + 1e-20)
            gradient = np.gradient(log_Fl_raw)
            peaks, _ = find_peaks(np.abs(gradient), height=0.3, distance=10)

            if len(peaks) == 0:
                spike_idx = len(Fl_raw)
                Fl_final = Fl_raw
            else:
                spike_idx = peaks[0]
                Fl_valid = Fl_raw[:spike_idx + 1]
                Fl_tail = extrapolate_tail(l_range, Fl_raw, spike_idx)
                Fl_final = np.concatenate([Fl_valid, Fl_tail[1:]])

            # Corrected sigma expressions
            sigma_raw = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * ((current_r / R)**l_range) * Fl_raw
            sigma = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * ((current_r / R)**l_range) * Fl_final

            ax.plot(l_range, sigma, color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')
            ax.plot(l_range, sigma_raw, color=colors[i], linestyle=':', alpha=0.6)  # dotted = original

        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-15, top=1e-4)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.tight_layout()
    plt.show()




def plot_sst_smoothened():
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

            # SMOOTHENED FILTER
            Fl_smooth = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2)) * (1 / np.sqrt(l_range))

            # Compute sigma using the smoothened filter
            sigma_smooth = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * (current_r / R)**l_range * Fl_smooth

            # Plot
            ax.plot(l_range, sigma_smooth, color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')

        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Smoothened Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-16, top=1e-4)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity SST measurement with LRI+ACC smoothened')
    plt.tight_layout()
    plt.show()


def plot_sst_v3():
    # Parameters
    h = np.array([100000, 150000, 200000, 250000, 300000])  # satellite altitudes (m)
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    d_fixed = 220000  # satellite separation (m)
    f = 1  # observation frequency (Hz) 1 measurement per second

    SOL_DURATION = 88775.244  # seconds in one Martian sol

    t_5_sols = 5 * SOL_DURATION
    t_30_sols = 30 * SOL_DURATION
    t_1_year = 668.6 * SOL_DURATION

    epsilon = 1e-8
    colors = ['b', 'g', 'r', 'purple', 'orange', 'cyan']
    l_range = np.arange(1, 401)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    durations = [
        (t_5_sols, "5 Sols", axs[0], 50),
        (t_30_sols, "30 Sols", axs[1], 90),
        (t_1_year, "1 Martian Year", axs[2], 90)
    ]

    for t, label, ax, l_max_compute in durations:
        N = f * t
        for i, height in enumerate(h):
            current_r = R + height
            current_n = n[i]

            gamma = 2 * np.arctan((d_fixed/2)/current_r)

            # Truncate or compute only up to max degree
            current_l_range = l_range[:l_max_compute]

            with np.errstate(divide='ignore', invalid='ignore'):
                Fl_raw = current_l_range / (np.sqrt(2 + 4 * current_l_range + 4 * current_l_range**2))
                Fl_raw *= (1 / np.abs(np.sin(current_l_range * gamma / 2) + 1e-10))

            # Compute sigma for the specific duration
            if label == "1 Martian Year":
                # Adjust sigma based on longer observation time and scaling
                sigma = (epsilon / np.sqrt(f * t_1_year)) * (1 / (current_n * current_r)) * (current_r / R)**current_l_range * Fl_raw
            else:
                # Compute sigma for 5-day and 30-day solutions
                sigma = (epsilon / np.sqrt(N)) * (1 / (current_n * current_r)) * (current_r / R)**current_l_range * Fl_raw

            # Prepare full sigma with exponential extrapolation
            full_sigma = np.zeros(len(l_range))
            full_sigma[:len(sigma)] = sigma
            
            # Exponential extrapolation for the rest of the range
            if len(sigma) > 1:
                # Compute the growth rate from the last few computed points
                last_segment = sigma[-5:]
                log_last_segment = np.log(last_segment)
                growth_rate = np.polyfit(current_l_range[-5:], log_last_segment, 1)[0]
                
                # Extrapolate with this growth rate
                extrapolation_start = len(sigma)
                extrapolation_x = l_range[extrapolation_start:] - l_range[extrapolation_start-1]
                full_sigma[extrapolation_start:] = sigma[-1] * np.exp(growth_rate * extrapolation_x)

            # Plot full line up to computed points
            ax.plot(l_range[:l_max_compute], full_sigma[:l_max_compute], 
                    color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')
            
            # Plot dashed extrapolation for the rest
            ax.plot(l_range[l_max_compute:], full_sigma[l_max_compute:], 
                    color=colors[i], linestyle='--', alpha=0.7)

        # Add RMS signal power law
        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l integrated ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-16, top=1e-4)
        ax.set_xlim(0, 400)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity SST measurement with LRI+ACC integrated')
    plt.tight_layout()
    plt.show()


def plot_sst_dist_satellites():
    import numpy as np
    import matplotlib.pyplot as plt

    # Mars constants
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)

    # Parameters
    height = 200000  # satellite altitude (m)
    d_values = np.array([50000, 100000, 150000, 200000, 250000])  # separations (m)
    f = 1  # observation frequency (Hz)
    epsilon = 1e-8
    l_range = np.arange(1, 401)
    l_max_compute = 50

    SOL_DURATION = 88775.244  # seconds in one Martian sol

    durations = [
        (5 *SOL_DURATION, "5 Sols"),
        (30 * SOL_DURATION, "30 Sols"),
        (668.6 * SOL_DURATION, "1 Martian Year")
    ]

    # Color palette
    colors = ['#1E90FF', '#32CD32', '#FF4500', '#9400D3', '#FF1493']

    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    for ax, (t_duration, label) in zip(axs, durations):
        for j, d_fixed in enumerate(d_values):
            current_r = R + height
            n = np.sqrt(mu / current_r**3)
            gamma = 2 * np.arctan((d_fixed / 2) / current_r)

            current_l_range = l_range[:l_max_compute]
            with np.errstate(divide='ignore', invalid='ignore'):
                Fl_raw = current_l_range / np.sqrt(2 + 4 * current_l_range + 4 * current_l_range**2)
                Fl_raw *= 1 / (np.abs(np.sin(current_l_range * gamma / 2)) + 1e-10)

            sigma = (epsilon / np.sqrt(f * t_duration)) * (1 / (n * current_r)) * (current_r / R)**current_l_range * Fl_raw

            # Prepare full sigma with extrapolation
            full_sigma = np.zeros(len(l_range))
            full_sigma[:len(sigma)] = sigma

            if len(sigma) > 1:
                last_segment = sigma[-5:]
                log_last_segment = np.log(last_segment)
                growth_rate = np.polyfit(current_l_range[-5:], log_last_segment, 1)[0]
                extrap_x = l_range[l_max_compute:] - current_l_range[-1]
                full_sigma[l_max_compute:] = sigma[-1] * np.exp(growth_rate * extrap_x)

            ax.plot(l_range[:l_max_compute], full_sigma[:l_max_compute],
                    color=colors[j], linestyle='-', linewidth=2,
                    label=f'{d_fixed / 1000:.0f} km sep')
            ax.plot(l_range[l_max_compute:], full_sigma[l_max_compute:],
                    color=colors[j], linestyle=':', linewidth=1.5, alpha=0.7)

        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})', fontsize=13)
        ax.set_xlabel('Spherical Harmonic Order (l)', fontsize=12)
        ax.set_yscale('log')
        ax.set_ylim(1e-16, 1e-4)
        ax.set_xlim(0, 400)
        ax.grid(True, which="both", ls="-", alpha=0.2)
        ax.legend(fontsize=9)

    axs[0].set_ylabel('Measurement Accuracy (σ)', fontsize=12)
    fig.suptitle('Sensitivity SST Measurement with LRI+ACC\nAltitude = 200 km, Varying Satellite Separation', fontsize=15)
    plt.tight_layout()
    plt.show()



def plot_sst_acc_vs_cai():
    # Constants
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    h = 200000  # satellite altitude (m)
    r = R + h
    n = np.sqrt(mu / r**3)

    d_fixed = 150000  # satellite separation (m)
    f = 1  # 1 Hz observation frequency

    # Time durations
    SOL_DURATION = 88775.244  # seconds in one Martian sol
    t_5_sols = 5 * SOL_DURATION
    t_30_sols = 30 * SOL_DURATION
    t_1_year = 668.6 * SOL_DURATION

    # Parameters to compare
    epsilon_vals = [1e-8, 1e-10]
    epsilon_labels = [r'LRI+ACC ($\epsilon = 10^{-8}$)', r'LRI+CAI ($\epsilon = 10^{-10}$)']
    epsilon_colors = ['blue', 'red']

    l_range = np.arange(1, 401)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    durations = [
        (t_5_sols, "5 Sols", axs[0], 50),
        (t_30_sols, "30 Sols", axs[1], 90),
        (t_1_year, "1 Martian Year", axs[2], 90)
    ]

    for t, label, ax, l_max_compute in durations:
        N = f * t
        gamma = 2 * np.arctan((d_fixed / 2) / r)
        current_l_range = l_range[:l_max_compute]

        with np.errstate(divide='ignore', invalid='ignore'):
            Fl_raw = current_l_range / (np.sqrt(2 + 4 * current_l_range + 4 * current_l_range**2))
            Fl_raw *= (1 / np.abs(np.sin(current_l_range * gamma / 2) + 1e-10))

        for eps, eps_label, eps_color in zip(epsilon_vals, epsilon_labels, epsilon_colors):
            sigma = (eps / np.sqrt(N)) * (1 / (n * r)) * (r / R)**current_l_range * Fl_raw

            # Fill full sigma with extrapolation
            full_sigma = np.zeros(len(l_range))
            full_sigma[:len(sigma)] = sigma

            if len(sigma) > 1:
                last_segment = sigma[-5:]
                log_last_segment = np.log(last_segment)
                growth_rate = np.polyfit(current_l_range[-5:], log_last_segment, 1)[0]
                extrapolation_start = len(sigma)
                extrapolation_x = l_range[extrapolation_start:] - l_range[extrapolation_start - 1]
                full_sigma[extrapolation_start:] = sigma[-1] * np.exp(growth_rate * extrapolation_x)

            ax.plot(l_range[:l_max_compute], full_sigma[:l_max_compute],
                    color=eps_color, linestyle='-', linewidth=2, label=eps_label)
            ax.plot(l_range[l_max_compute:], full_sigma[l_max_compute:],
                    color=eps_color, linestyle='--', alpha=0.6)

        # Plot observed signal power law
        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(1e-19, 1e-4)
        ax.set_xlim(0, 400)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity of SST Measurement (LRI+ACC) for Different ε')
    plt.tight_layout()
    plt.show()


# Run the function
plot_sst()
#plot_sst_v2()              #tried to integrate
#plot_sst_smoothened()      #here Fl is smoothened
plot_sst_v3()                #Integrated following Milano 
plot_sst_dist_satellites()
plot_sst_acc_vs_cai()