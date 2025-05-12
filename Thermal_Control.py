import numpy as np
from scipy.integrate import quad
import math 

# Operating Ranges:

# Payload:
# LRI: 27 to 30 C (operating)
# CAI: 10 C (operating, max grad. is 0.1 C)

# TT&C & CDHS
# Antennas: -100 to 100 C
# Antenna gimballs: -40 to 80 C
# C&DH Box baseplates: -20 to 60 C

# Structure
# Cylinder r = 1.5 m, l = 4.6 m

# ADCS
# Star Sensors: 0 to 30 C

# Propulsion
# Tank: 15 to 50

# Power - solar array area 8.4
# Solar array: 10 to 25 C

#standard variables
R_mars = 3390       #km
h_orbit = 200       #km
sig = 5.67 * 10**(-8)            #boltzman constant
T_Mars = 209.8
S_mars = 586.2
Albedo_mars = 0.25

#hot variables
alpha_s = 0.28                  #solar absorptivity, 127-mum silvered teflon after 5 years
eps = 0.75                      #emissivity 127-mum silvered teflon
Q_loss = 0                      #Heat lost or gained by surrouding structure
Q_d = 1150                      #internal energy dissipation
T_int = 30 + 273                    #Internal temperature
solar_incident = 0              #solar incidence angle

# #cold variables
# alpha_s = 0.08                  #solar absorptivity, 127-mum silvered teflon BOL
# eps = 0.79                      #emissivity, 127-mum silvered teflon
# Q_loss =                        #Heat lost or gained by surrouding structure
# Q_d = 950                       #internal energy dissipation


#Horizontal F factor as function of nu
def F12H(nu):
    # Define the integrand function
    def integrand(phi):
        A = 1 + (1 + nu)**2 - 2 * (1 + nu) * np.cos(phi) - np.sin(phi)**2
        B = 1 + (1 + nu)**2 - 2 * (1 + nu) * np.cos(phi)
        numerator = np.sqrt(A) * (np.sin(phi) * np.cos(phi) - np.sin(phi)**3)
        denominator = B**2
        return numerator / denominator

    # Upper limit of integration
    upper_limit = np.arccos(1 / (1 + nu))

    # Compute the integral
    result, _ = quad(integrand, 0, upper_limit)

    return 2 * result

#q_ir calculation
nu = h_orbit / R_mars
F = F12H(nu)
M_IR = sig * T_Mars**4
q_ir = F * M_IR #W/m^2

#q_sol
q_sol = S_mars * math.cos(solar_incident)

#q_A
q_A = S_mars * Albedo_mars * F

#Radiator area (hot data)
A_rad = (Q_d - Q_loss) / (eps * sig * T_int^4 - alpha_s * (q_sol + q_A) + eps * q_ir)

#Cold temperature (cold data)
sigT4 = (Q_d - Q_loss) / (eps * A_rad) + alpha_s / eps *(q_sol +q_A) + q_ir
T_cold = sigT4**(1/4) / sig

