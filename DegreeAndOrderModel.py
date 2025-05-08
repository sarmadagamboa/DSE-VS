import numpy as np
import matplotlib.pyplot as plt

def run_sim(h, delta, epsilon, N, setup):
    # Constants
    R = 3.3895e6 # Mars radius in meters
    mu = 4.282837e13 # m^3/s^2
    r = R + h
    n = np.sqrt(mu/(r**3))
    l = np.arange(1, 200, dtype=np.float64)
    if setup == 'SST':
        gamma = np.arctan(delta/r)
        F_SST = l/(np.sqrt(2+4*l+4*l**2))*(1/abs(np.sin(l*gamma/2)))
        Error_SST = F_SST*(epsilon/np.sqrt(N))*(1/(n*r))*(r/R)**l
        return l, Error_SST
    elif setup == 'single':
        F_single = l/(np.sqrt(1/2+l+l**2))
        Error_single = F_single*(epsilon/np.sqrt(N))*(1/(n*r))*(r/R)**l
        return l, Error_single
    elif setup == 'SGG':
        F_SGG = 1/(np.sqrt(2*(1+2*l+2*l**2)*(2+2*l+l**2)))
        Error_SGG = F_SGG*(epsilon/np.sqrt(N))*(1/(n**2))*(r/R)**l
        return l, Error_SGG

    
if __name__ == "__main__":
     # type = 'SST'
    setup = 'SST'
    h = 200e3 # Orbital altitude in meters
    delta = 100e3 # Separation between satellites in meters
    epsilon = 1e-5 # Measurement uncertainty in meters
    N = 27 # Number of observations
    l, error_single = run_sim(h, delta, epsilon, N, 'single')   
    l, error_SST = run_sim(h, delta, epsilon, N, 'SST')
    l, error_SGG = run_sim(h, delta, epsilon, N, 'SGG')

    # N = 50
    # l, error_single2 = run_sim(h, delta, epsilon, N, 'single')

    # N = 5
    # l, error_single3 = run_sim(h, delta, epsilon, N, 'single')

k = 1.1e-5 ##not completely sure about this value, but it is supposed to be the a constant based on the body i.e. Mars

RMS = k/l**2

plt.plot(l, error_SST, label='SST Error', color='blue')
plt.plot(l, error_single, label='Single Error at N = 500', color='red')
# plt.plot(l, error_single2, label='Single Error at N = 50', color='orange')
# plt.plot(l, error_single3, label='Single Error at N = 5', color='purple')
plt.plot(l, error_SGG, label='SGG Error', color='green')
plt.plot(l, RMS, label='Observed Signal Power Law', color='black', linestyle='--')
plt.xlabel('Spherical Harmonic Degree (l)')
plt.ylabel('Error (m)')
plt.title('Error in Spherical Harmonic Coefficients vs Degree')
plt.yscale('log')
#plt.xscale('log')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()