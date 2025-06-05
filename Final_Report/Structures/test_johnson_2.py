import math
import numpy as np
import matplotlib.pyplot as plt

E = 70e9
rho = 2800
height = 4.5  # meters
launch_load = 144814  # Newtons
min_weight = np.inf
best_config = None

def calculate_r(n_st, rho, height, t): 
    l0, l1, l2 = [0.03, 0.03, 0.03] # side flange, vertical web, top plate

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

# Try multiple thicknesses and number of stringers
thicknesses = np.linspace(0.002, 0.006, 9)  # m
n_stringers = range(8, 65, 4)

for t in thicknesses:
    for n_st in n_stringers:
        r_gyr, A_total, weight = calculate_r(n_st, rho, height, t)
        lambda_ratio = height / r_gyr
        sigma_euler = (np.pi**2 * E) / (lambda_ratio**2)
        total_euler_load = sigma_euler * A_total * n_st

        if total_euler_load >= launch_load and weight < min_weight:
            min_weight = weight
            best_config = {
                "thickness": t,
                "n_st": n_st,
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
