import math
import numpy as np
import matplotlib.pyplot as plt
from Structures2_optimal_current import Material, Stringer
from itertools import product

def calculate_r_hat(n_st, height, t, l0, l1, l2): 
     # side flange, vertical web, top plate

    # Segment areas
    a_side = t * l0
    a_vert = t * l1
    a_top = t * l2
    a_corner = t * t

    # y-positions of centroids from base
    y_side = t / 2
    y_vert = l1 / 2 + t
    y_top = l1 + t + t / 2
    y_corner_bottom = t / 2
    y_corner_top = t + l1 + t/2

    A_total = 2 * a_side + 2 * a_vert + a_top + 4 * a_corner
    y_bar = (2 * a_side * y_side + 2 * a_vert * y_vert + a_top * y_top + 
             2 * a_corner * y_corner_bottom + 2 * a_corner * y_corner_top) / A_total

    I_total = 0
    I_total += 2 * ((1/12) * l0 * t**3 + a_side * (y_bar - y_side)**2)
    I_total += 2 * ((1/12) * t * l1**3 + a_vert * (y_bar - y_vert)**2)
    I_total += (1/12) * l2 * t**3 + a_top * (y_bar - y_top)**2
    I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_top)**2)
    I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_bottom)**2)

    r_gyr = (I_total / A_total) ** 0.5
    weight = A_total * n_st * height * rho 
    return r_gyr, A_total, weight


def calculate_r_z(n_st, rho, height, t,l0, l1, l2): 
     # flange, web, flange

    # Segment areas
    a_side = t * l0     # bottom flange
    a_vert = t * l1     # web
    a_top = t * l2      # top flange
    a_corner = t * t

    # Assume web is vertical, bottom flange at y=0
    y_side = t / 2
    y_vert = l1 / 2 + t
    y_top = l1 + t + t / 2
    y_corner_top = l1 + t + t / 2
    y_corner_bottom = t/2

    A_total = a_side + a_vert + a_top + 2*a_corner
    y_bar = (a_side * y_side + a_vert * y_vert + a_top * y_top + a_corner * y_corner_top + a_corner * y_corner_bottom) / A_total

    # Moments of inertia about the centroid
    I_total = 0
    I_total += (1/12) * l0 * t**3 + a_side * (y_bar - y_side)**2
    I_total += (1/12) * t * l1**3 + a_vert * (y_bar - y_vert)**2
    I_total += (1/12) * l2 * t**3 + a_top * (y_bar - y_top)**2
    I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_top)**2)
    I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_bottom)**2)

    r_gyr = (I_total / A_total) ** 0.5
    weight = A_total * n_st * height * rho

    return r_gyr, A_total, weight


def determine_stringer_Euler(E, rho, height):
    # --- Constants ---
      # meters
    launch_load = 144814  # Newtons
    min_weight = np.inf
    best_config = None
    

    # Try multiple thicknesses and number of stringers
    thicknesses = np.linspace(0.002, 0.006, 90)  # m
    n_stringers = range(8, 65, 4)
    lengths = np.linspace(0.015, 0.05, 4)  # Try 15 mm to 50 mm for l0, l1, l2
    '''materials = [
        Material(name="Aluminium 7075", E=70e9, rho=2800, s_yld=448e6, s_ult=524e6, cost_kg=5.0, manuf=0.9, alpha=0.8, n=0.6, nu=0.334, alpha_therm = 23.6e-6, min_realistic_t = 0.0015),  # Aluminium 7075,
        Material(name="Titanium Ti-6Al-4V", E=116e9, rho=4420, s_yld=880e6, s_ult=950e6, cost_kg=30.0, manuf=0.7, alpha=0.75, n=0.15, nu=0.34, alpha_therm = 8.6e-6, min_realistic_t = 0.0020),  # Titanium Ti-6Al-4V
        Material(name="Stainless Steel", E=230e9, rho=7850, s_yld=250e6, s_ult=460e6, cost_kg=1.5, manuf=0.5, alpha=0.7, n=0.25, nu=0.3, alpha_therm = 17.2e-6, min_realistic_t = 0.0030),   # Stainless Steel
        Material(name="CFRP", E=150e9, rho=1600, s_yld=800e6, s_ult=1100e6, cost_kg=90.0, manuf=0.4, alpha=0.6, n=0.05, nu=0.2, alpha_therm = 0.5e-6, min_realistic_t = 0.0015)]  # CFRP
    stringer_types = [] 
    #stringer_types = [Stringer(type='hat', thickness=0.005, lengths=[0.03,0.03,0.03], material=material, manuf_stringer = 0.7),
            Stringer(type='Z', thickness=0.0028, lengths=[0.05,0.05,0.05], material=material, manuf_stringer = 0.9),
            Stringer(type='I', thickness=0.006, lengths=[0.03,0.03], material=material, manuf_stringer = 1),
            ] 
    '''

    for t in thicknesses:
        for n_st in n_stringers:
            for l0 in lengths:
                for l1 in lengths:
                    for l2 in lengths:
                        r_gyr, A_total, weight = calculate_r_z(n_st, rho, height, t, l0, l1, l2)
                        lambda_ratio = height / r_gyr
                        sigma_euler = (np.pi**2 * E) / (lambda_ratio**2)
                        total_euler_load = sigma_euler * A_total * n_st

                        if total_euler_load >= launch_load and weight < min_weight:
                            min_weight = weight
                            best_config = {
                                    #"Material": material.name,
                                    #"Stringer type": stringer.type,
                                    "Thickness of stringer": t,
                                    "n_st": n_st,
                                    "l0": l0,
                                    "l1": l1,
                                    "l2": l2,
                                    "r_gyr": r_gyr,
                                    "A_total": A_total,
                                    "Euler load": total_euler_load,
                                    "Weight": weight
                            }

    # Output result
    print("Best Configuration (Minimum Weight that passes Euler buckling):")
    for k, v in best_config.items():
        if isinstance(v, float):
            print(f"{k}: {v:.6f}")
        else:
            print(f"{k}: {v}")
    
    return best_config



def find_combinations(target_sum, num_digits):
    # Generate all combinations of 4 digits (from 0 to target_sum)
    
    configs = []
    if target_sum % 2 != 0:
        raise ValueError("Total number of stringers must be even to maintain symmetry.")
    
    half = target_sum // 2
    for s0 in range(half + 1):
        s1 = half - s0
        configs.append([s0, s1, s0, s1])
    return configs



def test_johnson_2f(): 
    ####### MODIFY INPUTS 
    E = 70e9
    rho = 2800
    height = 3

    best_config = determine_stringer_Euler(E, rho, height) #material and height fixed 
    target_sum = (best_config["n_st"])
    num_digits = 4
    configs = find_combinations(target_sum, num_digits)

    # Output the combinations
    #for config in configs:
        #print(config)
    
    return configs, height 
#determine_stringer_Euler()

configs, height = test_johnson_2f()
print(configs)
        