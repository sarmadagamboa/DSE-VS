import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm

def plot_comparative_measurements():
    # Parameters
    altitude = 200000  # fixed satellite altitude (m) - 200km
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + altitude
    n = np.sqrt(mu / r**3)
    
    # Satellite separations - now technique-specific
    separations = {
        'LRI+ACC (GRACE-FO)': 100000,  
        'LRI+CAI': 80000,            
        'Doppler': 0,             
        'QGG': 0                  
    }
        
    SOL_DURATION = 88775.244  # seconds in one Martian sol
    t_3yr_earth = 3 * 365 * 24 * 60 * 60  # 3 Earth years in seconds
    
    # Define measurement techniques with different epsilon values and markers
    measurement_techniques = [
        {
            'name': 'LRI+ACC (GRACE-FO)',
            'epsilon': 1e-8,
            'frequency': 1,
            'marker': 'o',
            'linestyle': '-',
            'markersize': 8
        },
        {
            'name': 'LRI+CAI',
            'epsilon': 1e-10,
            'frequency': 1,
            'marker': 's',
            'linestyle': '--',
            'markersize': 8
        },
        {
            'name': 'Doppler',
            'epsilon': 5e-6,
            'frequency': 0.0167,
            'marker': '^',
            'linestyle': '-.',
            'markersize': 8
        },
        {
            'name': 'QGG',
            'epsilon': 1e-10,
            'frequency': 1,
            'marker': 'x',
            'linestyle': ':',
            'markersize': 8
        }
    ]
    
    # Generate consistent line colors (all in red range)
    colors = ['darkred', 'red', 'orangered', 'darkorange']
    
    l_range = np.arange(1, 401, dtype=np.float64)
    
    plt.figure(figsize=(12, 8))
    
    # Plot each measurement type with different markers and colors
    for i, technique in enumerate(measurement_techniques):
        epsilon = technique['epsilon']
        f = technique['frequency']
        N = f * t_3yr_earth  # Sample count based on frequency
        
        if technique['name'] in ['LRI+ACC (GRACE-FO)', 'LRI+CAI']:
            # SST-hl formula with technique-specific separation
            d = separations[technique['name']]
            gamma = 2 * np.arctan((d / 2) / r)
            
            with np.errstate(divide='ignore', invalid='ignore'):
                Fl = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
                Fl *= (1 / np.abs(np.sin(l_range * gamma / 2)))
            
            sigma = (epsilon / np.sqrt(N)) * (1 / (n * r)) * ((r / R)**l_range) * Fl
            
        elif technique['name'] == 'QGG':
            # QGG formula with technique-specific separation
            d = separations[technique['name']]
            gamma = 2 * np.arctan((d / 2) / r)
            
            Fl = 1 / np.sqrt(2 * (1 + 2*l_range + 2*l_range**2) * (2 + 2*l_range + l_range**2))
            sigma = (epsilon / np.sqrt(N)) * (1 / n**2) * (r / R)**l_range * Fl
            
        elif technique['name'] == 'Doppler':
            # Doppler formula with technique-specific separation
            d = separations[technique['name']]
            gamma = 2 * np.arctan((d / 2) / r)
            
            Fl = l_range / np.sqrt(1/2 + l_range + l_range**2)
            sigma = (epsilon / np.sqrt(N)) * (1 / (n * r)) * (r / R)**l_range * Fl
        
        # Determine marker frequency based on line length
        marker_step = max(10, len(l_range) // 20)  # Show about 20 markers per line
        
        # Create mask for markers to avoid cluttering
        marker_mask = np.zeros_like(l_range, dtype=bool)
        marker_mask[::marker_step] = True
        
        # Plot line with line styles only (no markers on the line)
        plt.plot(
            l_range, 
            sigma, 
            color=colors[i], 
            linestyle=technique['linestyle'],
            linewidth=2,
            label=f"{technique['name']} (ε={epsilon:.0e}, f={f} Hz)"
        )
    
    # Plot the reference signal power law
    rms_signal = 8.5e-5 / l_range**2
    plt.plot(l_range, rms_signal, 'k--', linewidth=2, label='Observed Signal Power Law')
    
    # Format plot
    #plt.title(f'Measurement Sensitivity Comparison\n(Altitude: 200km, LRI+ACC: 80km, LRI+CAI: 100km, Duration: 3 years)')
    plt.xlabel('Spherical Harmonic Order (l)')
    plt.ylabel('Measurement Accuracy (σ)')
    plt.yscale('log')
    plt.ylim(bottom=1e-18, top=1e-4)
    plt.xlim(0, 400)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    
    # Create a custom legend with line styles, markers and colors
    legend = plt.legend(loc='best', fontsize=10, framealpha=0.9)
    legend.get_frame().set_edgecolor('black')
    
    # No colorbar needed anymore
    
    plt.tight_layout()
    plt.show()

# Run the function
plot_comparative_measurements()