import numpy as np

def calculate_sustainability():

    Li_ion_CO2 = 79.4                   # kg CO2-eq / kg bat
    Wh_kg = 150                         # LiCoO2
    Wh_l = 337                          # LiCoO2
    PV_CO2 = 412.9                      # kg CO2-eq / m2 arrays
    Xe_CO2 = 28.92                      # kg CO2-eq / kg Xe
    MMH_N2O4_CO2 = 0.96                 # kg CO2-eq / kg MMH (with N2O4)
    MMH_CO2 = 4.13*(10**3)*0.13         # kg CO2-eq / kg MMH  
    
    battery_sizes = np.array([0.0025, 0.0027, 0.003, 0.001]) # m^3, concept 1, 2, 3, 4
    battery_sizes *= 1000               # liters
    battery_Wh = battery_sizes * Wh_l
    battery_mass = battery_Wh / Wh_kg
    battery_CO2 = battery_mass * Li_ion_CO2

    array_areas = np.array([12, 15, 17.3, 9.4])
    array_CO2 = array_areas * PV_CO2

    Xe_mass = np.array([0, 8.8, 8.84, 0])
    MMH_N2O4_mass = np.array([343, 340.5, 342, 0])
    MMH_mass = np.array([0, 0, 0, 48.5])

    propellant_CO2 = Xe_mass * Xe_CO2 + MMH_N2O4_mass * MMH_N2O4_CO2 + MMH_mass * MMH_CO2

    rare_EM_use = np.array([2, 4, 3, 1])

    sustainability_impact = battery_CO2 + array_CO2 + propellant_CO2

    sustainability_ranks = sustainability_impact.argsort()
    sustainability_scores = np.zeros_like(sustainability_impact, dtype=int)
    sustainability_scores[sustainability_ranks] = [4, 3, 2, 1]

    rare_ranks = rare_EM_use.argsort()
    rare_scores = np.zeros_like(rare_EM_use, dtype=int)
    rare_scores[rare_ranks] = [4, 3, 2, 1]

    avg_scores = (sustainability_scores + rare_scores) / 2

    labels = ["LRI-ACC", "LRI-CAI", "QGG", "DT"]

    print("=== Sustainability Impact Summary ===")
    for i in range(len(labels)):
        print(f"{labels[i]}:")
        print(f"  → Total CO₂-eq impact:  {sustainability_impact[i]:.2f} kg")
        print(f"        → CO₂-eq (battery):    {battery_CO2[i]:.2f} kg")
        print(f"        → CO₂-eq (solar array): {array_CO2[i]:.2f} kg")
        print(f"        → CO₂-eq (propellant):  {propellant_CO2[i]:.2f} kg")
        print(f"  → Unsustainable material use:       {rare_EM_use[i]} units\n")
        print(f"  → Combined average score: {avg_scores[i]:.2f}")
        print(f"        → Sustainability score:   {sustainability_scores[i]}")
        print(f"        → Rare earth score:       {rare_scores[i]}")

    return sustainability_impact

calculate_sustainability()