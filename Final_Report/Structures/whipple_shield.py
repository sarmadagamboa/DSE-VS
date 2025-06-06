def whipple_rear_wall_thickness(
    m_p, velocity_kms, s, rho_p=2700, rho_b=2800,
):
    """
    https://ntrs.nasa.gov/api/citations/20140001399/downloads/20140001399.pdf
    Calculate required rear wall thickness in a Whipple shield setup.
    Units:
        - debris_d: particle diameter, m
        - velocity_kms: km/s
        - spacing: m - bumper to wall spacing 
        - densities: kg/m³
        - empirical constraints: Empirical constants (from NASA-HDBK-6003 for Al-Al systems):
    Returns:
        - Required rear wall thickness (m)
    """
    c = 0.055 * (rho_p*rho_b)**(1/6)
    t_r = c * m_p * velocity_kms / (s**(2))
    return t_r

m_p = 1e-2   
velocity_kms = 20    
s = 10          

t_rear = whipple_rear_wall_thickness(m_p, velocity_kms, s)
print(f"Required rear wall thickness: {t_rear*1e3:.2f} mm")

