import numpy as np
import matplotlib.pyplot as plt

def plot_qgg():
    # Parameters
    h = np.array([100000, 150000, 200000, 250000, 300000])  # satellite altitudes (m)
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    f = 0.5  # observation frequency (Hz) 1Hz max possible (use smaller)

    SOL_DURATION = 88775.244  # seconds in one Martian sol

    t_5_sols = 5 * SOL_DURATION
    t_30_sols = 30 * SOL_DURATION
    t_1_sols = 668.6 * SOL_DURATION

    epsilon = 1e-10 
    colors = ['b', 'g', 'r', 'purple', 'orange', 'cyan']
    l_range = np.arange(1, 401,dtype=np.float64)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    durations = [(t_5_sols, "5 Sols", axs[0],50), (t_30_sols, "30 Sols", axs[1],90), (t_1_sols, "1 Martian Year ", axs[2],90)]

    for t, label, ax, cutoff in durations:
        N = f * t
        for i, height in enumerate(h):
            current_r = R + height
            current_n = n[i]

            Fl = 1/ np.sqrt(2 * (1 + 2*l_range + 2*l_range**2 ) * (2 + 2*l_range + l_range**2))
            sigma = (epsilon/np.sqrt(N))*(1/current_n**2)*(current_r/R)**l_range*Fl
            
            ax.plot(l_range, sigma, color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')

        # Plot the reference signal power law
        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-12, top=1e-4)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity SST measurement with QGG')
    plt.tight_layout()
    plt.show()

def plot_doppler():
    # Parameters
    h = np.array([100000, 150000, 200000, 250000, 300000])  # satellite altitudes (m)
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + h
    n = np.sqrt(mu / r**3)

    f = 1  # observation frequency (Hz) ?? for QGG

    SOL_DURATION = 88775.244  # seconds in one Martian sol

    t_5_sols = 5 * SOL_DURATION
    t_30_sols = 30 * SOL_DURATION
    t_1_sols = 668.6 * SOL_DURATION

    epsilon = 1e-5 #Ka-Band
    colors = ['b', 'g', 'r', 'purple', 'orange', 'cyan']
    l_range = np.arange(1, 401,dtype=np.float64)

    fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    durations = [(t_5_sols, "5 Sols", axs[0],50), (t_30_sols, "30 Sols", axs[1],90), (t_1_sols, "1 Martian Year ", axs[2],90)]

    for t, label, ax, cutoff in durations:
        N = f * t
        for i, height in enumerate(h):
            current_r = R + height
            current_n = n[i]

            Fl = l_range/np.sqrt(1/2+ l_range + l_range**2)
            sigma = (epsilon/np.sqrt(N))*(1/(current_n*current_r))*(current_r/R)**l_range*Fl
            
            ax.plot(l_range, sigma, color=colors[i], linestyle='-', label=f'{height/1000:.0f} km')

        # Plot the reference signal power law
        rms_signal = 8.5e-5 / l_range**2
        ax.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')

        ax.set_title(f'Measurement Accuracy vs. l ({label})')
        ax.set_xlabel('Spherical Harmonic Order (l)')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e-12, top=1e-4)
        ax.grid(True)
        ax.legend()

    axs[0].set_ylabel('Measurement Accuracy (σ)')
    plt.suptitle('Sensitivity SST measurement with Doppler')
    plt.tight_layout()
    plt.show()

    

plot_qgg()
plot_doppler()
