from Final_Report.Structures.Structure2_optimal_current2 import structural_mass_wpanel as calc_struct_mass
import matplotlib.pyplot as plt

import math
g = 9.80665
"""
Calculates the structural nass 
    
    structural_mass_itercode = structural_mass_wpanel(sc_mass = 921.51, dim_height = 3, dim_length = 1.2, dim_width = 1.7)[0]

calc_prop_mass(dry_mass, Transfer_DeltaV, Onorbit_deltaV)
    returns the new propellant mass
"""

"""
calc_struct_mass(wet_mass, structural_setup)
    returns the structural mass
"""

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

    M_press = inputs["Propulsion setup"]["final_press"] * v_prop_biprop / (inputs["Propulsion setup"]['gas_const']* inputs["Propulsion setup"]['storage_press'])
    v_press = M_press * inputs["Propulsion setup"]['gas_const'] * inputs["Propulsion setup"]['storage_temp'] / inputs["Propulsion setup"]['storage_press']

    # Check tank capacities
    # Electric propellant check
    max_electric_prop = inputs["Propulsion setup"]["Electric_tank_volume"] * inputs["Propulsion setup"]["Electric_propellant_density"]
    electric_fits = (m_prop_stkeeping_electric <= max_electric_prop)

    # Fuel tank check
    max_fuel_capacity = inputs["Propulsion setup"]["Fuel_tank_volume"]
    fuel_fits = (v_fuel_required <= max_fuel_capacity)

    # Oxidizer tank check
    max_oxidizer_capacity = inputs["Propulsion setup"]["Oxidizer_tank_volume"]
    oxidizer_fits = (v_oxidizer_required <= max_oxidizer_capacity)

    # Pressurant tank check
    max_pressurant_capacity = inputs["Propulsion setup"]["Pressurant_tank_volume"]
    pressurant_fits = (v_press <= max_pressurant_capacity)

    all_tanks_fit = electric_fits and fuel_fits and oxidizer_fits and pressurant_fits
    total_prop_mass = fuel_mass+oxidiser_mass+m_prop_stkeeping_electric+M_press
    return total_prop_mass
    

def check_tank_capacities(dry_mass, inputs):
    """Check if the tanks have enough capacity for the propellant mass.
    """
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

    M_press = inputs["Propulsion setup"]["final_press"] * v_prop_biprop / (inputs["Propulsion setup"]['gas_const']* inputs["Propulsion setup"]['storage_press'])
    v_press = M_press * inputs["Propulsion setup"]['gas_const'] * inputs["Propulsion setup"]['storage_temp'] / inputs["Propulsion setup"]['storage_press']

    # Check tank capacities
    # Electric propellant check
    max_electric_prop = inputs["Propulsion setup"]["Electric_tank_volume"] * inputs["Propulsion setup"]["Electric_propellant_density"]
    electric_fits = (m_prop_stkeeping_electric <= max_electric_prop)

    # Fuel tank check
    max_fuel_capacity = inputs["Propulsion setup"]["Fuel_tank_volume"]
    fuel_fits = (v_fuel_required <= max_fuel_capacity)

    # Oxidizer tank check
    max_oxidizer_capacity = inputs["Propulsion setup"]["Oxidizer_tank_volume"]
    oxidizer_fits = (v_oxidizer_required <= max_oxidizer_capacity)

    # Pressurant tank check
    max_pressurant_capacity = inputs["Propulsion setup"]["Pressurant_tank_volume"]
    pressurant_fits = (v_press <= max_pressurant_capacity)

    all_tanks_fit = electric_fits and fuel_fits and oxidizer_fits and pressurant_fits
    tank_fits = [electric_fits, fuel_fits, oxidizer_fits, pressurant_fits]
    tank_vols = [max_electric_prop, max_fuel_capacity, max_oxidizer_capacity, max_pressurant_capacity]
    tank_caps = [m_prop_stkeeping_electric, v_fuel_required, v_oxidizer_required, v_press]
    return all_tanks_fit, tank_fits, tank_vols, tank_caps


def calculate_wet_mass(inputs, dry_mass_margin=1.1):
    print(inputs["mass"]["Propellant_mass"])
    return sum(value for key, value in inputs["mass"].items() if key != "Propellant_mass") * dry_mass_margin + inputs["mass"]["Propellant_mass"]


def calculate_dry_mass(inputs, dry_mass_margin=1.1):
    return sum(value for key, value in inputs["mass"].items() if key != "Propellant_mass") * dry_mass_margin


def iteration_loop(inputs, dry_mass_margin=1.1, values_close_percent=0.5):
    """Iterate until the wet mass converges to a stable value.
    """

    dim_height, dim_length, dim_width = inputs["Structural setup"]['dim_height'], inputs["Structural setup"]['dim_length'], inputs["Structural setup"]['dim_width']

    prop_mass_evolution = []
    wet_mass_evolution = []

    inputs["wet_mass_guess"] = calculate_dry_mass(inputs) + inputs["Structural_mass_guess"] * dry_mass_margin + inputs["Propellant_mass_guess"]
    print(inputs["wet_mass_guess"])

    wet_mass_evolution.append(inputs["wet_mass_guess"])
    prop_mass_evolution.append(inputs["Propellant_mass_guess"])

    inputs["mass"]["Structural_mass"] = calc_struct_mass(inputs["wet_mass_guess"], dim_height, dim_length, dim_width)[0]
    dry_mass = calculate_dry_mass(inputs)
    inputs["mass"]["Propellant_mass"] = calc_prop_mass(dry_mass, inputs)

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

        inputs["mass"]["Structural_mass"] = calc_struct_mass(calculate_wet_mass(inputs), dim_height, dim_length, dim_width)[0]
        dry_mass = calculate_dry_mass(inputs)
        inputs["mass"]["Propellant_mass"] = calc_prop_mass(dry_mass, inputs)

        wet_mass_evolution.append(calculate_wet_mass(inputs))
        prop_mass_evolution.append(inputs["mass"]["Propellant_mass"])

        if abs(wet_mass_evolution[-1] - wet_mass_evolution[-2]) / wet_mass_evolution[-2] * 100 < values_close_percent:
            values_close = True
        else:
            values_close = False
    
    inputs["tank_check"]["all_tanks_fit"], inputs["tank_check"]["tank_fits"], inputs["tank_check"]["tank_vols"], inputs["tank_check"]["tank_caps"] = check_tank_capacities(inputs["mass"]["Propellant_mass"], inputs)

    return wet_mass_evolution, prop_mass_evolution, inputs




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
        "Fuel_tank_volume" :0.165,  # m^3
        "Oxidizer_tank_volume" : 0.165,  # m^
        'Pressurant_tank_volume': 0.032,  # m^3

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


    # Assume wet_mass_evolution and prop_mass_evolution are already defined

    # Create a range for iteration numbers
    iterations = range(len(wet_mass_evolution))

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(iterations, wet_mass_evolution, label='Wet Mass', marker='o')
    plt.plot(iterations, prop_mass_evolution, label='Propellant Mass', marker='x')

    # Add labels and title
    plt.xlabel('Iteration Number')
    plt.ylabel('Mass (kg)')
    plt.title('Evolution of Wet Mass and Propellant Mass')
    plt.legend()
    plt.grid(True)

    # Display the plot
    plt.tight_layout()
    plt.show()
