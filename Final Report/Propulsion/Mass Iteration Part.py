
import math
g = 9.80665

if __name__ == "__main__":
    inputs = {}

    inputs["mass"] = {
        "Payload_mass": 142,  # kg
        "ADCS_mass": 36.7,  # kg
        "TTC_mass": 75,  # kg
        "CDHS_mass": 10,  # kg
        "Thermal_mass": 29.45,  # kg
        "Power_mass": 41.36,  # kg
        "Propulsion_dry_mass": 134,  # kg
    }

    inputs["Propulsion setup"] = {
        "Electric_isp":1500,  # s
        "Electric_prop_margin": 0.20,
        "Electric_propellant_density": 1350,  # kg/m^3 (Xenon)
        "Electric_tank_volume": 2 * 0.005,

        "Biprop_isp": 321, # s
        "Biprop_prop_margin": 0.15,
        "Fuel_density" :880,  # kg/m^3 (MMH)
        "Oxidizer_density": 1370,  # kg/m^3 (MON-3)
        "Ox_fuel_ratio" :1.65,
        "Fuel_tank_volume" :0.198,  # m^3
        "Oxidizer_tank_volume" : 0.198,  # m^
        'Pressurant_tank_volume': 0.040, #0.032 m^3

        'final_press': 2000000,  # Pa
        'storage_temp': 300,  # K
        'storage_press': 20000000,  # Pa
        'gas_const': 2077,  # J/kg*K
    }

    inputs["Structural setup"] = {
        'dim_height': 3,  # m
        'dim_length': 1.2,  # m
        'dim_width': 1.7,  # m
        }

    inputs["Propellant_mass_guess"] = 0 # kg
    inputs["Structural_mass_guess"] = 86.48 # kg
    inputs["Insertion_DeltaV"] = 1308.14  # m/s -- total deltaV of the Mars insertion
    inputs["Onorbit_DeltaV"] = 56+196.21  # m/s -- total deltaV of the spacecraft after insertion
    inputs["Capture_Time"] = 45 #minutes
    inputs["Stationkeeping_Time"] = 20 * 427 + 365 * 24 * 60 / 10 #minutes

def calc_prop_mass(dry_mass, inputs):
    m_beforesk = dry_mass * math.exp(inputs["Onorbit_DeltaV"] / (inputs["Propulsion setup"]["Electric_isp"] * g))
    m_prop_electric_nomargin = m_beforesk - dry_mass
    m_prop_stkeeping_electric = m_prop_electric_nomargin * (1 + inputs["Propulsion setup"]["Electric_prop_margin"])

    m_beforecapture = m_beforesk * math.exp(inputs["Insertion_DeltaV"] / (inputs["Propulsion setup"]["Biprop_isp"] * g))
    m_prop_capture_nomargin = m_beforecapture - m_beforesk
    m_prop_capture_biprop =  m_prop_capture_nomargin*(1+inputs["Propulsion setup"]["Biprop_prop_margin"])

    fuel_mass = m_prop_capture_biprop / (1 + inputs["Propulsion setup"]["Ox_fuel_ratio"])
    oxidiser_mass = inputs["Propulsion setup"]["Ox_fuel_ratio"] * fuel_mass

    v_fuel_required = fuel_mass/inputs["Propulsion setup"]["Fuel_density"]
    v_oxidizer_required = oxidiser_mass /inputs["Propulsion setup"]["Oxidizer_density"]
    v_prop_biprop = v_fuel_required + v_oxidizer_required

    M_press = inputs["Propulsion setup"]["final_press"] * v_prop_biprop / (inputs["Propulsion setup"]['gas_const']* inputs["Propulsion setup"]['storage_temp'])
    v_press = M_press * inputs["Propulsion setup"]['gas_const'] * inputs["Propulsion setup"]['storage_temp'] / inputs["Propulsion setup"]['storage_press']

    # Check tank capacities
    # Electric propellant check
    max_electric_prop = inputs["Propulsion setup"]["Electric_tank_volume"] * inputs["Propulsion setup"]["Electric_propellant_density"]
    electric_fits = m_prop_stkeeping_electric <= max_electric_prop

    # Fuel tank check
    max_fuel_capacity = inputs["Propulsion setup"]["Fuel_tank_volume"]
    fuel_fits = v_fuel_required <= max_fuel_capacity

    # Oxidizer tank check
    max_oxidizer_capacity = inputs["Propulsion setup"]["Oxidizer_tank_volume"]
    oxidizer_fits = v_oxidizer_required <= max_oxidizer_capacity

    # Pressurant tank check
    max_pressurant_capacity = inputs["Propulsion setup"]["Pressurant_tank_volume"]
    pressurant_fits = (v_press <= max_pressurant_capacity)

    all_tanks_fit = electric_fits and fuel_fits and oxidizer_fits and pressurant_fits

    if all_tanks_fit:
        # Return total propellant mass if everything fits
        total_prop_mass = fuel_mass+oxidiser_mass+m_prop_stkeeping_electric+M_press
        return print(total_prop_mass)
    else:
        # Return None if tanks don't fit
        return print(max_oxidizer_capacity, v_oxidizer_required, max_fuel_capacity, v_fuel_required, v_press, max_pressurant_capacity)

dry_mass = 680

print(calc_prop_mass(dry_mass,inputs))