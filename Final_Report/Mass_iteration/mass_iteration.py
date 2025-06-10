from Final_Report.Structures.Structure2_optimal_current2 import structural_mass_wpanel

from Final_Report. import ... as calc_prop_mass ###THIS IS BULLSHIT
import math
g = 9.80665
"""
Calculates the structural nass 
    
    structural_mass = structural_mass_wpanel(sc_mass = 1063)

calc_prop_mass(dry_mass, Transfer_DeltaV, Onorbit_deltaV)
    returns the new propellant mass
"""


def calculate_wet_mass(inputs, dry_mass_margin=1.1):

    m_beforesk = dry_mass * math.exp(inputs["Onorbit_DeltaV"] / (inputs["Electric_isp"]["Propulsion setup"] * g))
    m_prop_electric_nomargin = m_beforesk - dry_mass
    m_prop_stkeeping_electric = m_prop_electric_nomargin * (1 + inputs["Electric_prop_margin"]["Propulsion setup"])

    m_beforecapture = m_beforesk * math.exp(inputs["Insertion_DeltaV"] / (inputs["Biprop_isp"]["Propulsion setup"] * g))
    m_prop_capture_nomargin = m_beforecapture - m_beforesk
    m_prop_capture_biprop =  m_prop_capture_nomargin*(1+inputs["Biprop_prop_margin"]["Propulsion setup"])

    fuel_mass = m_prop_capture_biprop / (1 + inputs["Ox_fuel_ratio"]["Propulsion setup"])
    oxidiser_mass = inputs["Ox_fuel_ratio"]["Propulsion setup"] * fuel_mass

    v_fuel_required = fuel_mass/inputs["Fuel_density"]["Propulsion setup"]
    v_oxidizer_required = oxidiser_mass /inputs["Oxidizer_density"]["Propulsion setup"]
    v_prop_biprop = v_fuel_required + v_oxidizer_required

    M_press = inputs["final_press"]["Propulsion setup"] * v_prop_biprop / (inputs['gas_const']["Propulsion setup"]* inputs['storage_press']["Propulsion setup"])
    v_press = M_press * inputs['gas_const']["Propulsion setup"] * inputs['storage_temp']["Propulsion setup"] / inputs['storage_press']["Propulsion setup"]

    # Check tank capacities
    # Electric propellant check
    max_electric_prop = inputs["Electric_tank_volume"]["Propulsion setup"] * inputs["Electric_propellant_density"]["Propulsion setup"]
    electric_fits = m_prop_stkeeping_electric <= max_electric_prop

    # Fuel tank check
    max_fuel_capacity = inputs["Fuel_tank_volume"]["Propulsion setup"]
    fuel_fits = v_fuel_required <= max_fuel_capacity

    # Oxidizer tank check
    max_oxidizer_capacity = inputs["Oxidizer_tank_volume"]["Propulsion setup"]
    oxidizer_fits = v_oxidizer_required <= max_oxidizer_capacity

    # Pressurant tank check
    max_pressurant_capacity = inputs["Pressurant_tank_volume"]["Propulsion setup"]
    pressurant_fits = (v_press <= max_pressurant_capacity)

    all_tanks_fit = electric_fits and fuel_fits and oxidizer_fits and pressurant_fits

    if all_tanks_fit:
        # Return total propellant mass if everything fits
        total_prop_mass = fuel_mass+oxidiser_mass+m_prop_stkeeping_electric+M_press
        return total_prop_mass
    else:
        # Return None if tanks don't fit
        return None


    return sum(value for key, value in inputs["mass"].items() if key != "Propellant_mass") * dry_mass_margin + inputs["mass"]["Propellant_mass"]


def calculate_dry_mass(inputs, dry_mass_margin=1.1):
    return sum(value for key, value in inputs["mass"].items() if key != "Propellant_mass") * dry_mass_margin


def iteration_loop(inputs, dry_mass_margin=1.1, values_close_percent=0.5):
    """Iterate until the wet mass converges to a stable value.
    """
    prop_mass_evolution = []
    wet_mass_evolution = []

    inputs["wet_mass_guess"] = calculate_dry_mass(inputs) + inputs["Structural_mass_guess"] * dry_mass_margin + inputs["Propellant_mass_guess"]
    print(inputs["wet_mass_guess"])

    wet_mass_evolution.append(inputs["wet_mass_guess"])
    prop_mass_evolution.append(inputs["propellant_mass_guess"])

    inputs["mass"]["Structural_mass"] = calc_struct_mass(inputs["wet_mass_guess"], inputs["Structural setup"])
    dry_mass = calculate_dry_mass(inputs)
    inputs["mass"]["Propellant_mass"] = calc_prop_mass(dry_mass, inputs["Insertion_DeltaV"], inputs["Onorbit_DeltaV"])

    wet_mass_evolution.append(calculate_wet_mass(inputs))
    prop_mass_evolution.append(inputs["mass"]["Propellant_mass"])

    if abs(wet_mass_evolution[-1] - wet_mass_evolution[-2]) / wet_mass_evolution[-2] * 100 < values_close_percent:
        values_close = True
    else:
        values_close = False

    counter = 0

    while not values_close:
        counter += 1
        if counter > 100:
            print("Warning: Mass iteration loop exceeded 100 iterations. Check inputs or convergence criteria.")
            break

        inputs["mass"]["Structural_mass"] = calc_struct_mass(calculate_wet_mass(inputs), inputs["Structural setup"])
        dry_mass = calculate_dry_mass(inputs)
        inputs["mass"]["Propellant_mass"] = calc_prop_mass(dry_mass, inputs["Insertion_DeltaV"], inputs["Onorbit_DeltaV"])

        wet_mass_evolution.append(calculate_wet_mass(inputs))
        prop_mass_evolution.append(inputs["mass"]["Propellant_mass"])

        if abs(wet_mass_evolution[-1] - wet_mass_evolution[-2]) / wet_mass_evolution[-2] * 100 < values_close_percent:
            values_close = True
        else:
            values_close = False
    
    return wet_mass_evolution, prop_mass_evolution, inputs




if __name__ == "__main__":
    inputs = {}

    inputs["mass"] = {
        "Payload_mass": 102,  # kg
        "ADCS_mass": 65,  # kg
        "TTC_mass": 49,  # kg
        "CDHS_mass": 10,  # kg
        "Thermal_mass": 29,  # kg
        "Power_mass": 108,  # kg
        "Propulsion_dry_mass": 134,  # kg
    }

    inputs["Propulsion setup"] = {
        "Electric_isp":1500,  # s
        "Electric_prop_margin": 0.20,
        "Electric_propellant_density": 1350,  # kg/m^3 (Xenon)
        "Electric_tank_volume": 2 * 0.004,

        "Biprop_isp": 321, # s
        "Biprop_prop_margin": 0.15,
        "Fuel_density" :880,  # kg/m^3 (MMH)
        "Oxidizer_density": 1370,  # kg/m^3 (MON-3)
        "Ox_fuel_ratio" :1.65,
        "Fuel_tank_volume" :0.165,  # m^3
        "Oxidizer_tank_volume" : 0.165,  # m^
        'Pressurant_tank_volume': 0.032,  # m^3

        'final_press': 2000000,  # Pa
        'storage_temp': 300,  # K
        'storage_press': 20000000,  # Pa
        'gas_const': 2077,  # J/kg*K
    }

    inputs["Structural setup"] = {
        "REVIEW": 0
        }

    inputs["Propellant_mass_guess"] = 0 # kg
    inputs["Structural_mass_guess"] = 86.48 # kg
    inputs["Insertion_DeltaV"] = 1249  # m/s -- total deltaV of the Mars insertion
    inputs["Onorbit_DeltaV"] = 0.125 * 290.1 + 196.21  # m/s -- total deltaV of the spacecraft after insertion
    inputs["Capture_Time"] = 45 #minutes
    inputs["Stationkeeping_Time"] = 20 * 290.1 + 365 * 24 * 60 / 10 #minutes

#add valuesclose_percent to inputs
#add dry_mass_margin to inputs

    wet_mass_evolution, prop_mass_evolution, outputs = iteration_loop(inputs)

    final_wet_mass = wet_mass_evolution[-1]
    final_prop_mass = prop_mass_evolution[-1]
    final_struct_mass = outputs["mass"]["Structural_mass"]
 
    print(f"Final wet mass: {final_wet_mass} kg")
    print(f"Final propellant mass: {final_prop_mass} kg")
    print(f"Final structural mass: {final_struct_mass} kg")