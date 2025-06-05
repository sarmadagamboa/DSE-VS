import math
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

def calculate_r(): 
    t = 0.0005
    l0, l1, l2 = [0.01, 0.01, 0.01] # side flange, vertical web, top plate

    # Segment areas
    a_side = t * l0
    a_vert = t * l1
    a_top = t * l2
    a_corner = t * t  # corner reinforcement patches

    # y-positions of centroids from base (assuming flat base)
    y_side = t / 2                     # base flange center
    y_vert = l1 / 2 + t               # web center
    y_top = l1 + t + t / 2            # top plate center
    y_corner_bottom = t / 2              # corners (if placed just above base)
    y_corner_top = t + l1 + t/2

    # Area-weighted centroid (neutral axis)
    A_total = 2 * a_side + 2 * a_vert + a_top + 4 * a_corner
    y_bar = (2 * a_side * y_side + 2 * a_vert * y_vert + a_top * y_top + 2 * a_corner * y_corner_bottom + 2 * a_corner * y_corner_top) / A_total

    # got here - Moments of inertia about the centroidal axis
    I_total = 0
    # Sides (base flanges)
    I_total += 2 * ((1/12) * l0 * t**3 + a_side * (y_bar - y_side)**2)
    # Verticals
    I_total += 2 * ((1/12) * t * l1**3 + a_vert * (y_bar - y_vert)**2)
    # Top plate
    I_total += (1/12) * l2 * t**3 + a_top * (y_bar - y_top)**2
    # Corners
    I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_top)**2)
    I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_bottom)**2)

    # Radius of gyration
    r_gyr = (I_total / A_total) ** 0.5
    return r_gyr
  
#hat stringer crippling stress: 
stringer_stress = 221878525.32905626
height = 4.5
E = 70e9
yield_stress = 448e6
r_gyr = calculate_r()

lambda_vals = np.linspace(10, 1000, 1500)

#current_slenderness_ratio = height/r_gyr #to use in loop 
lambda_cr = math.sqrt(2 * math.pi**2 * E /stringer_stress) 

#For actual loading calculation 
# Johnson stress (only valid where λ < λ_cr)
sigma_johnson = np.where(
    lambda_vals < lambda_cr,
    stringer_stress*(1 - (stringer_stress * (lambda_vals)**2)/(4 * math.pi**2 * E) ),
    np.nan
)

# Euler stress (only valid where λ >= λ_cr)
sigma_euler = np.where(
    lambda_vals >= lambda_cr,
    np.pi**2 * E / (lambda_vals**2),
    np.nan
)

#for plotting
sigma_johnson_plot = stringer_stress*(1 - (stringer_stress * (lambda_vals)**2)/(4 * math.pi**2 * E) )
sigma_johnson_plot[sigma_johnson_plot < 0] = np.nan
sigma_euler_plot = np.pi**2 * E / (lambda_vals**2)
print(np.pi**2 * E / ((height/r_gyr)**2)) 

print(sigma_johnson_plot)
# Plotting
plt.figure(figsize=(10, 6))
plt.plot(lambda_vals, sigma_johnson_plot / 1e6, label="Johnson's crippling stress", lw=2)
plt.plot(lambda_vals, sigma_euler_plot / 1e6, label="Euler buckling stress", lw=2)
plt.axvline(x=lambda_cr, color='grey', linestyle='--', label=f'Critical λ ≈ {lambda_cr:.1f}')

plt.xlabel("Slenderness Ratio (λ = L/r)", fontsize=12)
plt.ylabel("Stress [MPa]", fontsize=12)
plt.title("Johnson vs Euler Buckling Stress vs Slenderness Ratio", fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


""" 
#s_johnson_crippling = stringer_stress(1 - (stringer_stress * (height*r_gyr)**2)/(4 * math.pi()**2 * E) )
#euler_buckling_yielding = math.pi()**2 * E * (r_gyr/height)**2
critical_slenderness_ratio = math.sqrt(2 * math.pi**2 * E /stringer_stress) 
current_slenderness_ratio = height/r_gyr #to use in loop 
#then based on that calculate limit 
print(critical_slenderness_ratio, current_slenderness_ratio)
#since we are above critical, then we have Euler buckling more likely, must be checked against yield strength 

#choose curve depending on scenario 
if current_slenderness_ratio > critical_slenderness_ratio: 
    print("Euler buckling critical. Check against yielding.")
    #euler curve 
    #check against yielding 
"""