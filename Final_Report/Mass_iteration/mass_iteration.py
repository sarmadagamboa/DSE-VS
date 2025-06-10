from Final_Report.Structures.Structure2_optimal_current2 import *

from Final_Report. import ... as calc_prop_mass ###THIS IS BULLSHIT

"""
Calculates the structural nass 
    
    structural_mass_itercode = structural_mass_wpanel(sc_mass = 921.51, dim_height = 3, dim_length = 1.2, dim_width = 1.7)[0]

calc_prop_mass(dry_mass, Transfer_DeltaV, Onorbit_deltaV)
    returns the new propellant mass
"""


def calculate_wet_mass(inputs, dry_mass_margin=1.1):
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
        "Propulsion_dry_mass": 72,  # kg
    }

    inputs["Structural setup"] = {
        "REVIEW": 0
        }

    inputs["Propellant_mass_guess"] = 0 # kg
    inputs["Structural_mass_guess"] = 86.48 # kg
    inputs["Insertion_DeltaV"] = 0  # m/s -- total deltaV of the Mars insertion
    inputs["Onorbit_DeltaV"] = 0  # m/s -- total deltaV of the spacecraft after insertion

#add valuesclose_percent to inputs
#add dry_mass_margin to inputs

    wet_mass_evolution, prop_mass_evolution, outputs = iteration_loop(inputs)

    final_wet_mass = wet_mass_evolution[-1]
    final_prop_mass = prop_mass_evolution[-1]
    final_struct_mass = outputs["mass"]["Structural_mass"]
 
    print(f"Final wet mass: {final_wet_mass} kg")
    print(f"Final propellant mass: {final_prop_mass} kg")
    print(f"Final structural mass: {final_struct_mass} kg")