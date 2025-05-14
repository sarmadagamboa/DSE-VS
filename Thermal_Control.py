import numpy as np
import math 

# ─── To do ─────────────────────────────────────────────
# automate temperature range selection
# add heaters
# add Q_LOSS

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

#T_COLD = max(T_COLD_LIST)  # Cold-case temperature [°C]
T_HOT = min(T_HOT_LIST)    # Hot-case temperature [°C]

# ─── Heaters ────────────────────────────────────────────────────────
# Patch heaters with solid-state thermostat have accuracy of less than 0.1 C 
A_ACC = 26 * 26 * 6               # CAI area [cm²]
A_LRI = 40 * 40 * 2 + 40 * 20 * 4 # LRI area [cm²]
A_PAYLOAD = A_ACC + A_LRI         # Payload area [cm²]
P_HEATER = 0.5  
U = 2.5                   # Heater power [W/cm²]
TCS_POWER = U * A_PAYLOAD/10000 * (28.5 - 19) # Power needed to heat the payload [W]

# ─── Physical Constants ────────────────────────────────────────────────────────
SIGMA = 5.67e-8        # Stefan–Boltzmann constant [W/m²·K⁴]

# ─── Mars & Orbit Parameters ─────────────────────────────────────────────────
R_MARS = 3389.5        # Mars radius [km]
H_ORBIT = 200.0        # Orbital altitude [km]
NU = H_ORBIT / R_MARS  # Dimensionless altitude ratio (a/R)
r = 1.5                # Cylinder radius [m]
l = 4.6                # Cylinder length [m]
Q_DISS = 800.0         # Internal dissipation [W]
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

def cylinder_area(r, l):
    """
    Total surface area of a cylinder with radius r and length l.
    """
    return 2 * math.pi * r * (l + r)

def compute_TCS_mass(A_rad, A_tot):
    """
    Returns the mass of the TCS [kg] based on its area [m²]:
    m = A_rad * ρ * t
    where ρ is the density of the radiator material and t is its thickness.
    """
    # Material properties (example values)
    rho_teflon = 2700.0  # Density of aluminum [kg/m³]
    t = 127e-6     # Thickness of the radiator [m]
    A_MLI = A_tot - A_rad
    rho_MLI = 0.73 # Density of MLI [kg/m²] for 15 layers
    return A_rad * rho_teflon * t + A_MLI * rho_MLI

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


def main():
    q_sol, q_alb, q_ir = compute_fluxes()
    A_rad = compute_radiator_area(q_sol, q_alb, q_ir)
    T_cold = compute_cold_temperature(A_rad, q_sol, q_alb, q_ir)
    A_tot = cylinder_area(r, l)  # Cylinder dimensions [m]
    TCS_mass = compute_TCS_mass(A_rad, A_tot)

    print(f"Dimensionless altitude ν = a/Rₘₐᵣₛ: {NU:.4f}")
    print(f"Radiator area       : {A_rad:.2f} m²")
    print(f"Cold-case temp      : {T_cold:.2f} °C")
    print(f"Fluxes [W/m²]       : q_sol={q_sol:.2f}, q_alb={q_alb:.2f}, q_ir={q_ir:.2f}")
    print(f"Total TCS mass      : {TCS_mass:.2f} kg")
    print(f"TCS Power           : {TCS_POWER:.2f} W")




if __name__ == "__main__":
    main()