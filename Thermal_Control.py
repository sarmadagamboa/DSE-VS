import numpy as np

# ─── To do ─────────────────────────────────────────────
# automate temperature range selection
# add heaters
# add Q_LOSS
# look into battery
# look into payload
# add MLI
# loook into transfer orbits

# ─── Mars Thermal Control System ─────────────────────────────────────────────
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

# ─── Operating temperature ranges ────────────────────────────────────────────────────────
T_PAYLOAD = [27.0, 30.0]       # Payload operating range [°C]
T_ADCS = [0.0, 30.0]           # ADCS operating range [°C]
T_TTC = [-20.0, 60.0]          # TT&C operating range [°C]
T_STRUCTURE = [-100.0, 100.0]  # Structure operating range [°C]
T_PROPULSION = [15.0, 50.0]    # Propellant tank operating range [°C]
T_POWER = [0.0, 20.0]          # Battery operating range [°C]

T_COLD_LIST = [T_PAYLOAD[0], T_ADCS[0], T_TTC[0], T_STRUCTURE[0], T_PROPULSION[0], T_POWER[0]]
T_HOT_LIST = [T_PAYLOAD[1], T_ADCS[1], T_TTC[1], T_STRUCTURE[1], T_PROPULSION[1], T_POWER[1]]
SUBSYSTEM_NAMES = ['Payload', 'ADCS', 'TT&C', 'Structure', 'Propulsion', 'Power']
T_HOT = 30.0 # min(T_HOT_LIST) # Hot-case temperature [°C]

# ─── Payload heater ────────────────────────────────────────────────────────
# Patch heaters with solid-state thermostat have accuracy of less than 0.1 C 
A_CAI = 0.52 * 0.52 * 6                           # CAI area [m²]
A_ACC = 0.26 * 0.26 * 6                           # ACC area [m²]
A_LRI = 0.40 * 0.40 * 2 + 0.40 * 0.20 * 4         # LRI area [m²]
A_QGG = 0.5 * 0.5 * 2 + 0.5 * 1  * 4              # QGG area [m²]     
A_PAYLOAD = A_ACC                         # Payload area [m²]
U = 2.5                                           # Heater power [W/m²]
T_PAYLOAD_AVG = (T_PAYLOAD[0] + T_PAYLOAD[1]) / 2 # Average payload temperature [°C]

# ─── Physical Constants ────────────────────────────────────────────────────────
SIGMA = 5.67e-8        # Stefan–Boltzmann constant [W/m²·K⁴]

# ─── Mars & Orbit Parameters ─────────────────────────────────────────────────
R_MARS = 3389.5        # Mars radius [km]
H_ORBIT = 200.0        # Orbital altitude [km]
NU = H_ORBIT / R_MARS  # Dimensionless altitude ratio (a/R)
w = 3.0                # Satellite width [m]
l = 3.0                # Satellite length [m]
h = 2.0                # Satellite height [m]
Q_DISS = 598.0         # Internal dissipation [W]
T_MARS = 209.8         # Mars effective temperature [K]
S_MARS = 586.2         # Solar constant at Mars [W/m²]
ALBEDO_MARS = 0.25     # Mars Bond albedo

# ─── Hot-Case Thermal Properties ─────────────────────────────────────────────
ALPHA_HOT = 0.28            # Solar absorptivity (127 μm Teflon, EOL)
EPSILON_HOT = 0.75          # Emissivity (127 μm Teflon)
Q_DISS_HOT = Q_DISS * 1.1   # Internal dissipation [W] with margin
Q_LOSS_HOT = 0.0            # Structural heat exchange [W]
T_INTERNAL = T_HOT + 273.15 # Internal set-point [K]
SOLAR_ANGLE = 0.0           # Incidence angle (rad)

# ─── Cold-Case Thermal Properties ────────────────────────────────────────────
ALPHA_COLD = 0.08          # Solar absorptivity (127 μm Teflon, BOL)
EPSILON_COLD = 0.79        # Emissivity (127 μm Teflon)
Q_DISS_COLD = Q_DISS * 0.9 # Internal dissipation [W] with margin
Q_LOSS_COLD = 0.0          # Structural heat exchange [W]
T_MIN_COLD = 15.0          # Minimum operating temperature (T_PROPULSION) [°C]

# ─── Functions ────────────────────────────────────────────────────────────────
def filter_subsystems_above_Tcold(T_cold, T_cold_list, names):
    """
    Returns a list of (subsystem_name, min_temp) tuples for those subsystems
    whose minimum operating temperature is higher than T_cold.
    """
    return [
        (name, t_min)
        for name, t_min in zip(names, T_cold_list)
        if t_min > T_cold
    ]

def compute_payload_heater_power(T_cold):
    """
    Returns the power needed to heat the payload [W] based on its area [m²]:
    P = U * A * (T_hot - T_cold)
    where U is the heater power [W/m²] and T_hot and T_cold are the hot and cold temperatures.
    """
    return U * A_PAYLOAD * (T_PAYLOAD_AVG - T_cold)

def compute_satellite_area(w, l, h):
    """
    Total surface area of a cylinder with radius r and length l.
    """
    return 2*(w*l+w*h+h*l)

def compute_TCS_mass(A_rad, A_tot):
    """
    Returns the mass of the TCS [kg] based on its area [m²]:
    m = A_rad * ρ * t
    where ρ is the density of the radiator material and t is its thickness.
    """
    # Material properties (example values)
    rho_teflon = 0.27  # Density of aluminum [kg/m²]
    A_MLI = A_tot - A_rad
    rho_MLI = 0.73 # Density of MLI [kg/m²] for 15 layers
    return A_rad * rho_teflon + A_MLI * rho_MLI

def compute_fluxes(view_factor=0.9):
    """
    Returns the three solar/planetary fluxes [W/m²] incident on a horizontal surface:
      - q_sol: direct solar
      - q_alb: Mars-albedo reflected
      - q_ir: Mars infrared
    """
    q_sol = S_MARS * np.cos(SOLAR_ANGLE)
    q_alb = S_MARS * ALBEDO_MARS * view_factor
    q_ir  = view_factor * SIGMA * T_MARS**4
    return q_sol, q_alb, q_ir


def compute_radiator_area(q_sol, q_alb, q_ir):
    """
    Radiator area [m²] needed in the hot case to reject Q_DISS_HOT:
    A = (Q_DISS_HOT - Q_LOSS_HOT) / [ ε·σ·T_int⁴ 
                                     - α·(q_sol+q_alb) 
                                     + ε·q_ir ]
    """
    numerator = Q_DISS_HOT - Q_LOSS_HOT
    denominator = (
        EPSILON_HOT * SIGMA * T_INTERNAL**4
        - ALPHA_HOT * (q_sol + q_alb)
        + EPSILON_HOT * q_ir
    )
    return numerator / denominator


def compute_cold_temperature(A_rad, q_sol, q_alb, q_ir):
    """
    Equilibrium cold-case temperature [°C] of that same radiator area:
    Solve σ·T⁴ = (Q_DISS_COLD - Q_LOSS_COLD)/(ε·A)
                 + (α/ε)·(q_sol+q_alb) + q_ir
    """
    term = (
        (Q_DISS_COLD - Q_LOSS_COLD) / (EPSILON_COLD * A_rad)
        + (ALPHA_COLD / EPSILON_COLD) * (q_sol + q_alb)
        + q_ir
    )
    T_cold_K = (term / SIGMA) ** 0.25
    return T_cold_K - 273.15

def compute_heater_power_coldcase(q_sol, q_alb, q_ir, A_rad):
    """
    Returns Q_heater [W] required in the cold case:
      Q_heater = Q_loss
               + A_r*( ε·σ·T_min^4   − [ α_s*(q_sol+q_alb) + ε·q_ir ] )
               − Q_d
    """
    # convert target temp to Kelvin
    T_min_K = T_MIN_COLD + 273.15

    # radiative power out at T_min
    rad_term = EPSILON_COLD * SIGMA * T_min_K**4

    # environmental loading
    env_term = ALPHA_COLD * (q_sol + q_alb) + EPSILON_COLD * q_ir

    # heater power
    P = Q_LOSS_COLD + A_rad * (rad_term - env_term) - Q_DISS_COLD
    if P < 0:
        print("T_cold > T_min, no heating power required.")
        return 0.0
    else:
        return P

def main():
    q_sol, q_alb, q_ir = compute_fluxes()
    A_tot = compute_satellite_area(w, l, h)
    A_rad = compute_radiator_area(q_sol, q_alb, q_ir)
    A_MLI = A_tot - A_rad
    T_cold = compute_cold_temperature(A_rad, q_sol, q_alb, q_ir)  
    TCS_mass = compute_TCS_mass(A_rad, A_tot)
    subsys_above = filter_subsystems_above_Tcold(T_cold, T_COLD_LIST, SUBSYSTEM_NAMES)
    P_PAYLOAD = compute_payload_heater_power(T_cold)
    P_GENERAL = compute_heater_power_coldcase(q_sol, q_alb, q_ir, A_rad)
    P_TOT = P_PAYLOAD + P_GENERAL

    # Print results
    print(f"Radiator area       : {A_rad:.2f} m²")
    print(f"MLI area            : {A_MLI:.2f} m²")
    print(f"Payload area        : {A_PAYLOAD:.2f} m²")
    print(f"Cold-case temp      : {T_cold:.2f} °C")
    print(f"Fluxes [W/m²]       : q_sol={q_sol:.2f}, q_alb={q_alb:.2f}, q_ir={q_ir:.2f}")
    print(f"Total TCS mass      : {TCS_mass:.2f} kg")
    print(f"Payload heater power: {P_PAYLOAD:.2f} W")
    print(f"General heater power: {P_GENERAL:.2f} W")
    print(f"Total power         : {P_TOT:.2f} W")
    # print(f"Max T allowed       : {T_HOT:.2f} °C")
    print("Subsystems requiring heating:")
    for name, t_min in subsys_above:
        print(f" - {name}: {t_min:.1f} °C > {T_cold:.2f} °C")

if __name__ == "__main__":
    main()