import numpy as np
import matplotlib.pyplot as plt
from marsdensity import plot_mars_density



def compute_atmospheric_drag(mars_altitude, spacecraft_velocity, spacecraft_mass, drag_coefficient, cross_sectional_area, delta_t, atmospheric_density):
    """
    Calculate delta V caused by atmospheric drag at a given altitude and density.

    spacecraft_velocity in km/s, converted inside.
    """
    velocity_m_s = spacecraft_velocity # km/s to m/s
    # print(atmospheric_density)
    drag_force = 0.5 * atmospheric_density* velocity_m_s**2 * drag_coefficient * cross_sectional_area
    #print(drag_force)
    deceleration = drag_force / spacecraft_mass  # m/s²
    delta_v_drag = deceleration * delta_t  # m/s over delta_t seconds

    
    return delta_v_drag


def plot_delta_v_vs_altitude(alt = 200, year=1, spacecraft_mass=850, drag_coefficient=2.6, cross_sectional_area=1.54, delta_t=2.2*365*24*3600):
    # Get heights and max densities from your function
    R_mars = 3390e3
    mu = 4.282837e13
    results = plot_mars_density(year)
    heights = results['height']
    densities = results['max_density']

    idx = np.abs(heights - alt).argmin()
    closest_alt = heights[idx]
    rho = densities[idx]
    # print(rho)
    if np.isnan(rho):
        raise ValueError(f"No atmospheric density data available for altitude {alt} km.")

    r = R_mars +closest_alt*1000
    v_orbit = np.sqrt(mu / r)
    # print(v_orbit)
    results = plot_mars_density(year)
    heights = results['height']
    densities = results['max_density']
    # print(f"THIS IS THE DENSITY: {densities}")
    delta_vs = []
    for h, rho in zip(heights, densities):
        if np.isnan(rho):
            delta_vs.append(np.nan)
        else:
            delta_v = compute_atmospheric_drag(
                mars_altitude=h,
                spacecraft_velocity=v_orbit,
                spacecraft_mass=spacecraft_mass,
                drag_coefficient=drag_coefficient,
                cross_sectional_area=cross_sectional_area,
                delta_t=delta_t,
                atmospheric_density=rho
            )
            delta_vs.append(delta_v)

    delta_vs = np.array(delta_vs)
    # print(f"Delta V values: {delta_vs}")
    plt.figure(figsize=(8, 6))
    plt.plot(delta_vs, heights, marker='o')
    plt.xlabel('Delta V required due to drag (m/s)')
    plt.ylabel('Altitude (km)')
    plt.title(f'Delta V vs Altitude (TES Year {year})')
    plt.grid(True)
    plt.gca().invert_yaxis()  # Optional: altitude decreases downward
    plt.show()


def orbital_decay(alt, year= 2, spacecraft_mass=700, drag_coefficient=2.6, cross_sectional_area=1.54):
    """
    Calculate the orbital decay due to atmospheric drag over a specified time period.
    """

    R_mars = 3390e3
    mu = 4.282837e13  # Mars gravitational parameter in m^3/s^2

    results = plot_mars_density(year)
    heights = results['height']
    densities = results['max_density']

    idx = np.abs(heights - alt).argmin()
    closest_alt = heights[idx]
    rho = densities[idx]
    print(rho)
    if np.isnan(rho):
        raise ValueError(f"No atmospheric density data available for altitude {alt} km.")

    r = R_mars +closest_alt*1000
    v_orbit = np.sqrt(mu / r) 
    print(v_orbit) # Orbital velocity in m/s
    drag_force = 0.5 * rho * v_orbit**2 * drag_coefficient * cross_sectional_area
    print(drag_force)
    s = 2*np.pi * r  # Circumference of the orbit
    delta_E = -drag_force * s 
    E = -spacecraft_mass * mu / (2*r)  # Initial orbital energy
    E_final = E + delta_E
    r_final = -spacecraft_mass * mu / (2 * E_final)
    #print(r_final)  # New radius after decay
    delta_r = (-r+r_final)/1000
    return delta_r 








# Example usage:
altitude = 212.48
plot_delta_v_vs_altitude(
    alt = altitude,
    year=2,# km/s example velocity
    spacecraft_mass=654.731,      # kg
    drag_coefficient=2.6,
    cross_sectional_area=1.45,  # m²
    delta_t=1.88*365*24*3600       # seconds duration
)


decay = orbital_decay(
    alt =altitude,
    year=2,
    spacecraft_mass=654.731,
    drag_coefficient=2.6,
    cross_sectional_area=1.45
)
print(f"Orbital decay at altitude {altitude} km: {decay} km")