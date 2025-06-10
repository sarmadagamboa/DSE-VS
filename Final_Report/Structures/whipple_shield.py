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

def whipple_mass(dim_width, dim_length): 
    m_p = 1e-2   
    velocity_kms = 20    
    s = 10          

    t_rear = whipple_rear_wall_thickness(m_p, velocity_kms, s)
    t_bumper = t_rear/5
    t_total = t_bumper + t_rear
    print(f"Required rear wall thickness: {t_rear*1e3:.2f} mm")
    print(f"Required bumper wall thickness: {t_bumper*1e3:.2f} mm")

    mass_whipple_shield = t_bumper * dim_width * dim_length * 2800
    print(mass_whipple_shield)
    return mass_whipple_shield
    

#since sufficient 1.5 mm wall monolithic wall sufficient - aluminium can act as the rear wall in a Whipple configuration,
#but does require bumper in front 
