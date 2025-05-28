from dataclasses import dataclass
import numpy as np

@dataclass
class Mission:
    """Mission parameters for solar‑array sizing at Mars (Mars LEO example)."""

    # Sub‑system loads (W)
    power_subsystems = {
        "ADCS":       212,
        "Structures":   0,
        "TT&C":       195,
        "CDHS":        20,
        "Power":        0,  # overwritten a few lines below
        "Thermal":     55,
        "Propulsion":   21,
        "Payload":     194,
    }
    mass_subsystems = {
        "ADCS":       62,
        "Structures":  149.7,
        "TT&C":       40,
        "CDHS":        10,
        "Power":       0,
        "Thermal":     29,
        "Propulsion":   56,
        "Payload":     32,
    }

    # power subsystem mass and power estimations based on literature
    power_req_subsystems: float = sum(power_subsystems.values())  # W
    power_subsystems["Power"] = 0.2 * power_req_subsystems       # 15margin for the power system itself

    # Power requirements
    power_req_eol: float = sum(power_subsystems.values())          # W at end‑of‑life
    power_margin: float = 0.0  #system margin                         


    # Orbit geometry
    period: float = 6529.2                  # s, orbital period
    eclipse: float = 2611.68                # s, eclipse duration
    sunlight: float = period - eclipse      # s, sunlight duration
    theta: float = 30                        # deg, theta angle (worst‑case)

    # Panel & pointing
    offpoint: float = 5.0          # deg off‑pointing for sun‑tracking arrays
    cell_temp: float = -40.0       # C worst‑case cell temperature
    temp_coeff: float = -0.25/100  # %/C power drop
    cell_eff: float = 0.30         # GaAs BOL efficiency
    packing_factor: float = 0.90   # geometric packing efficiency
    X_e: float = 0.65              # battery charge efficiency
    X_d: float = 0.85              # distribution efficiency

    # Solar environment
    solar_flux: float = 586.2      # W/m2
    P_o: float = cell_eff * solar_flux # GaAs, W/m2

    # Lifetime & degradation
    lifetime: float = 4.5          # years
    inherent_degradation: float = 0.85
    cell_degradation: float = 0.99 # annual degradation factor

    # Array configuration: 0 = body‑fixed, 1 = rotating / sun‑tracking
    array_type: int = 0

def kelly_cos(theta_deg):
    """Kelly cosine approximation for body fixed panels (>60 deg)."""

    KELLY = {0: 1.00, 30: 0.866, 50: 0.635, 60: 0.450, 80: 0.100, 85: 0.000}
    keys = np.array(sorted(KELLY))
    vals = np.array([KELLY[k] for k in keys])
    return np.interp(theta_deg, keys, vals)

def solar_array_sizing(m: Mission):
    """Return area, power_eol_density, power_bol, power_sa, incidence_factor."""

    # EOL power
    power_eol = m.power_req_eol * (1 + m.power_margin)

    # Required power during sunlight
    power_sa = ((power_eol * m.eclipse / m.X_e) + (power_eol * m.sunlight / m.X_d)) / m.sunlight

    # Incidence factor
    if m.array_type == 0:
        incidence_factor = kelly_cos(abs(m.theta))
    else:
        incidence_factor = kelly_cos(abs(m.offpoint))

    # BOL & EOL power density
    power_bol = m.P_o * m.inherent_degradation * incidence_factor
    L_d = m.cell_degradation ** m.lifetime
    power_eol_density = power_bol * L_d

    # Required array area
    area = (power_sa / power_eol_density)

    return area, power_eol, power_bol, power_sa, incidence_factor


def battery_sizing(m: Mission):
    """Compute battery parameters.
    """
    power_eol = m.power_req_eol * (1 + m.power_margin)
    DOD = 0.85
    efficiency_loss = 0.85
    E_bat = power_eol * (m.eclipse/3600) / (DOD * efficiency_loss)
    # Battery mass
    specific_power = 300  # W/kg, power book Li-ion
    energy_density = 400  # Wh/L, power book, Li-ion
    bat_mass = E_bat / specific_power
    bat_volume = E_bat / energy_density
    return bat_mass, bat_volume

if __name__ == "__main__":
    m = Mission()

    area, p_eol, p_bol, p_sa, inc = solar_array_sizing(m)
    bat_mass, bat_volume = battery_sizing(m)
    power_mass = (area * 2.06 + bat_mass + 0.071 * p_eol + 0.15)

    print("\n— Mars Solar-Array Quick Sizing —")
    print(f"Configuration               : {m.array_type}")
    print(f"Power required at EOL (W)   : {p_eol:,.2f}")
    print(f"BOL power density (W/m²)    : {p_bol:,.2f}")
    print(f"Subsystem power (W)         : {m.power_req_subsystems:,.2f}")
    print(f"Power while in sunlight (W) : {p_sa:,.2f}")
    print(f"Incidence factor            : {inc:.3f}")
    print(f"Array area (m²)             : {area:,.2f}")
    print(f"Power-system margin (W)     : {m.power_subsystems['Power']:,.2f}")
    print(f"Battery mass (kg)           : {bat_mass:,.2f}")
    print(f"Battery volume (L)          : {bat_volume:,.2f}")


    H = "\033[1;36m"   # bold-cyan
    R = "\033[0m"      # reset
    print(f"\n{H}===ITERATION PARAMETERS==={R}")
    print(f"{H}Power-system required power  : {m.power_subsystems['Power']:,.2f} W{R}")
    print(f"{H}Total power subsystem mass: {power_mass:,.2f} kg{R}")
    print(f"{H}Solar-array area     : {area:,.2f} m²{R}")
    print(f"{H}Total s/c power (W)            : {p_eol:,.2f}{R}\n")