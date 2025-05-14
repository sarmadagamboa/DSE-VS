import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import matplotlib.cm as cm

def plot_comparative_sst_measurements():
    # Parameters
    altitude = 200000  # fixed satellite altitude (m) - 200km
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + altitude
    n = np.sqrt(mu / r**3)
    
    d_fixed = 150000  # satellite separation (m) - 150km
        
    SOL_DURATION = 88775.244  # seconds in one Martian sol
    t_5_sols = 5 * SOL_DURATION  # Using 5 Sols for this comparison
    
    # Using only the original epsilon values (no multipliers)
    
    # Define measurement techniques with different epsilon values
    measurement_techniques = {
        'LRI+ACC (GRACE-FO)': {
            'epsilon': 1e-8,
            'frequency': 1,
            'color': 'blue'
        },
        'LRI+CAI': {
            'epsilon': 1e-10,
            'frequency': 0.5,
            'color': 'green'
        },
        'Doppler': {
            'epsilon': 1e-5,  #1e-5 K-band 
            'frequency': 0.0167,
            'color': 'red'
        },
        'QGG': {
            'epsilon': 1e-10,
            'frequency': 0.5,
            'color': 'purple'
        }
    }
    
    
    l_range = np.arange(1, 401, dtype=np.float64)
    gamma = 2 * np.arctan((d_fixed / 2) / r)
    
    plt.figure(figsize=(12, 8))
    
    # Plot each measurement type with different epsilon values
    for technique, params in measurement_techniques.items():
        epsilon = params['epsilon']
        f = params['frequency']
        color = params['color']
        N = f * t_5_sols  # Sample count based on frequency
       
        
        if technique in ['LRI+ACC (GRACE-FO)', 'LRI+CAI']:
            # SST-hl formula
            with np.errstate(divide='ignore', invalid='ignore'):
                Fl = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
                Fl *= (1 / np.abs(np.sin(l_range * gamma / 2)))
            
            sigma = (epsilon / np.sqrt(N)) * (1 / (n * r)) * ((r / R)**l_range) * Fl
            
        elif technique == 'QGG':
            # QGG formula
            Fl = 1 / np.sqrt(2 * (1 + 2*l_range + 2*l_range**2) * (2 + 2*l_range + l_range**2))
            sigma = (epsilon / np.sqrt(N)) * (1 / n**2) * (r / R)**l_range * Fl
            
        elif technique == 'Doppler':
            # Doppler formula
            Fl = l_range / np.sqrt(1/2 + l_range + l_range**2)
            sigma = (epsilon / np.sqrt(N)) * (1 / (n * r)) * (r / R)**l_range * Fl
        
        label = f'{technique} (ε={epsilon:.0e}, f={f} Hz)'
        
        plt.plot(l_range, sigma, color=color, linewidth=2, label=label)
    
    # Plot the reference signal power law
    rms_signal = 8.5e-5 / l_range**2
    plt.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')
    
    # Format plot
    plt.title(f'Measurement Sensitivity Comparison\n(Altitude: 200km, Separation: 150km, Duration: 5 Sols)')
    plt.xlabel('Spherical Harmonic Order (l)')
    plt.ylabel('Measurement Accuracy (σ)')
    plt.yscale('log')
    plt.ylim(bottom=1e-18, top=1e-4)
    plt.xlim(0, 400)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend(loc='best', fontsize=8)
    
    plt.tight_layout()
    plt.show()

# Run the function
plot_comparative_sst_measurements()