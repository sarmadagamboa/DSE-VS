#q_sun * alpha_MLI * A_solar + (q_ir + q_alb) * alpha_MLI * A_nad - Q_diss = eps_rad * sigma * T^4 * A_rad + eps_MLI * sigma * T^4 * A_MLI

# Constants
sigma = 5.67e-8  # Stefan-Boltzmann constant in W/(m²·K⁴)
A_nad = 1
A_solar = 1
A_MLI = 1
h_orbit = 212 # Height of orbit in km
R_Mars = 3390 # Mars radius in km
nu = h_orbit / R_Mars  # Ratio of orbit height to Mars radius
view_factor = 0.9 # View factor for at h = 212 km

# Hot values
alpha_MLI_hot = 0.0035
eps_MLI_hot = 0.013
eps_rad_hot = 0.75  # Emissivity of the radiator
Q_diss_hot = 1000
T_hot = 30 + 273.15 # Internal temperature in °C 
q_sun_hot = 715 # Solar flux in W/m² at periphilion
q_ir_hot = 300
q_alb_hot = 150

# Cold values
alpha_MLI_cold = 0.0023
eps_MLI_cold = 0.019
eps_rad_cold = 0.79  # Emissivity of the radiator
Q_diss_cold = 500 
q_sun_cold = 492 # Solar flux in W/m² at aphelion
q_ir_cold = 100
q_alb_cold = 50

def compute_radiator_area():
    """
    Radiator area [m²] needed in the hot case to reject Q_DISS:
    A_rad = (Q_diss - q_sun * alpha_MLI * A_solar - (q_ir + q_alb) * alpha_MLI * A_nad) / (eps_rad * sigma * T_hot^4)
    """
    numerator = q_sun_hot * alpha_MLI_hot * A_solar + (q_ir_hot + q_alb_hot) * alpha_MLI_hot * A_nad - Q_diss_hot - eps_MLI_hot * sigma * T_hot**4 * A_MLI
    denominator = eps_rad_hot * sigma * T_hot**4
    return numerator / denominator

if __name__ == "__main__":
    A_rad_hot = compute_radiator_area()
    print(f"Radiator area in hot case: {A_rad_hot:.2f} m²")
