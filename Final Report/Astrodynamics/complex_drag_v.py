import numpy as np
import matplotlib.pyplot as plt
from marsdensity import plot_mars_density

def compute_atmospheric_drag(mars_altitude, spacecraft_velocity, spacecraft_mass, drag_coefficient, cross_sectional_area, delta_t, atmospheric_density):
    """
    Calculate delta V caused by atmospheric drag at a given altitude and density.

    spacecraft_velocity in km/s, converted inside.
    """
    velocity_m_s = spacecraft_velocity * 1000  # km/s to m/s
    drag_force = 0.5 * atmospheric_density * velocity_m_s**2 * drag_coefficient * cross_sectional_area
    deceleration = drag_force / spacecraft_mass  # m/s²
    delta_v_drag = deceleration * delta_t  # m/s over delta_t seconds
    return delta_v_drag


def plot_delta_v_vs_altitude(year=1, spacecraft_velocity=4.7, spacecraft_mass=850, drag_coefficient=2.6, cross_sectional_area=1.54, delta_t=3.3*365*24*3600):
    # Get heights and max densities from your function
    results = plot_mars_density(year)
    heights = results['height']
    densities = results['max_density']

    delta_vs = []
    for h, rho in zip(heights, densities):
        if np.isnan(rho):
            delta_vs.append(np.nan)
        else:
            delta_v = compute_atmospheric_drag(
                mars_altitude=h,
                spacecraft_velocity=spacecraft_velocity,
                spacecraft_mass=spacecraft_mass,
                drag_coefficient=drag_coefficient,
                cross_sectional_area=cross_sectional_area,
                delta_t=delta_t,
                atmospheric_density=rho
            )
            delta_vs.append(delta_v)

    delta_vs = np.array(delta_vs)

    plt.figure(figsize=(8, 6))
    plt.plot(delta_vs, heights, marker='o')
    plt.xlabel('Delta V required due to drag (m/s)')
    plt.ylabel('Altitude (km)')
    plt.title(f'Delta V vs Altitude (TES Year {year})')
    plt.grid(True)
    plt.gca().invert_yaxis()  # Optional: altitude decreases downward
    plt.show()

# Example usage:
plot_delta_v_vs_altitude(
    year=1,
    spacecraft_velocity=4.7,  # km/s example velocity
    spacecraft_mass=850,      # kg
    drag_coefficient=2.6,
    cross_sectional_area=1.45,  # m²
    delta_t=3.3*365*24*3600                # seconds duration
)