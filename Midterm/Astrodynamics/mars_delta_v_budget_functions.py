
import numpy as np
from astropy.constants import G, M_earth, R_earth
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# Constants
M_mars = 6.4171e23 * u.kg
R_mars = 3389.5 * u.km
d_MS = 228e6  # km
d_ES = 149e6  # km
mu_earth = (G * M_earth).to(u.km**3 / u.s**2).value
mu_mars = (G * M_mars).to(u.km**3 / u.s**2).value

r_earth = R_earth.to(u.km).value
r_mars = R_mars.to(u.km).value

mu_sun = 1.32712440018e11  # km^3/s^2
rho_0 = 0.020  # Atmospheric density at the surface of Mars (kg/m³)
H = 11.1  # Scale height of the Martian atmosphere (km)
# Mission parameters
altitude_data = np.array([0, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400])
density_data = np.array([0.020, 0.009, 0.004, 0.002, 0.001, 0.0005, 0.0002, 5e-5, 1e-5, 1e-7, 1e-9, 1e-11, 1e-13])

# Interpolation function
density_interp = interp1d(altitude_data, density_data, kind='linear', fill_value='extrapolate')


def compute_launch_to_leo(leo_radius):
    return np.sqrt(mu_earth / leo_radius)

# def compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo):
#     transfer_a = (leo_radius + (d_MS + d_ES) + mars_orbit_radius) / 2
#     v_earth_orbit = np.sqrt(mu_sun / d_ES)
#     v_hel_LEO = np.sqrt(mu_sun * (2 / d_ES - 1 / transfer_a))
#     v_inf_leo = v_hel_LEO - v_earth_orbit
#     v_per_LEO = np.sqrt((2 * mu_earth / leo_radius) + v_inf_leo**2)
#     delta_v = np.abs(v_per_LEO - v_leo)
#     eccentricity = 1 - ((leo_radius + d_ES) / transfer_a)
#     return delta_v, transfer_a, eccentricity, v_inf_leo
def compute_dir_tranfer_injection(r_earth, mars_orbit_radius):
    transfer_a = ((d_MS + d_ES) + mars_orbit_radius)/ 2
    v_earth_orbit = np.sqrt(mu_sun / d_ES)
    v_hel_LEO = np.sqrt(mu_sun * (2 / d_ES - 1 / transfer_a))
    v_inf_leo = v_hel_LEO - v_earth_orbit
    v_per_LEO = np.sqrt((2 * mu_earth / r_earth) + v_inf_leo**2)
    delta_v = np.abs(v_per_LEO) 
    eccentricity = 1 - ((r_earth + d_ES) / transfer_a)
    return delta_v, transfer_a, eccentricity, v_inf_leo

def compute_capture_orbit(mars_orbit_radius):
    e_mars = 0.8537
    a_mars = mars_orbit_radius / (1 - e_mars)
    print(a_mars)
    r_apo_mars = 2 * a_mars - mars_orbit_radius
    v_mars_apo = np.sqrt(mu_mars * (-1 / a_mars + 2 / r_apo_mars))
    v_mars_per = np.sqrt(mu_mars * (-1 / a_mars + 2 / mars_orbit_radius))
    return r_apo_mars, v_mars_apo, v_mars_per, a_mars

def compute_mars_inclination_change(v_orbit):
    
    delta_v_incl = np.sqrt(2 * v_orbit**2 * (1 - np.cos(np.radians(93))))
    return delta_v_incl

def compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ):
    v_apo_mars = np.sqrt((2 * mu_mars / mars_orbit_radius) + v_inf_mars**2)
    return np.abs(v_mars_circ - v_apo_mars)
def compute_mars_period(mars_orbit_radius):
    T_mars = 2 * np.pi * np.sqrt(((mars_orbit_radius)**3) / mu_mars)
    return T_mars

def compute_station_keeping(years=4.5, delta_v_per_year=60):
    return delta_v_per_year * years / 1000  # km/s

def compute_deorbit():
    return 60 / 1000  # km/s


def compute_atmospheric_drag(mars_altitude, spacecraft_velocity, spacecraft_mass, drag_coefficient, cross_sectional_area, delta_t):
    """
    Compute the effect of atmospheric drag at a specified Martian altitude.

    Parameters:
        mars_altitude (float): Altitude above the Martian surface (km).
        spacecraft_velocity (float): Velocity of the spacecraft relative to the atmosphere (m/s).
        spacecraft_mass (float): Mass of the spacecraft (kg).
        drag_coefficient (float): Drag coefficient (dimensionless).
        cross_sectional_area (float): Cross-sectional area of the spacecraft (m²).

    Returns:
        delta_v_drag (float): Velocity change due to atmospheric drag (m/s).
    """
    h = mars_altitude * 1000  # Convert km to m

    # Compute atmospheric density at the given altitude
    # rho = rho_0 * np.exp(-h / (H * 1000))  # Density in kg/m³
    # T = -23.4-0.00222*h  # Temperature in K
    # p = 0.699 * np.exp(-0.00009*h)  # Pressure in Pa
    # rho = p/(0.1921*T)
    rho = 10**-11
    print(f"Atmospheric density at {mars_altitude} km: {rho} kg/m³")
    # Compute drag force
    drag_force = 0.5 * rho * (spacecraft_velocity*1000)**2 * drag_coefficient * cross_sectional_area  # Force in N
    print(f"Drag force: {drag_force} N")
    # Compute deceleration due to drag
    deceleration = drag_force / spacecraft_mass  # Acceleration in m/s²
    delta_v_drag = (deceleration * delta_t)/1000  # Velocity change in m/s
    # Output results
    # print(f"At altitude {mars_altitude} km:")
    # print(f"  Atmospheric density: {rho:.6f} kg/m³")
    # print(f"  Drag force: {drag_force:.3f} N")
    # print(f"  Deceleration: {deceleration:.6f} m/s²")
    # print(f"  Δv due to drag: {delta_v_drag:.3f} m/s")

    return delta_v_drag

def main(use_aerobraking, inclination_midcourse, leo_alt, mars_orbit_alt, spacecraft_mass, Cd, cross_sectional_area, delta_t):
    leo_radius = r_earth + leo_alt
    mars_orbit_radius = r_mars + mars_orbit_alt
    v_leo = np.sqrt(mu_earth / leo_radius)
    v_mars_circ = np.sqrt(mu_mars / mars_orbit_radius)
    # === CALCULATIONS ===
    # delta_v_launch = compute_launch_to_leo(leo_radius)
    #delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)
    delta_v_1, transfer_a, e, v_inf_leo = compute_dir_tranfer_injection(leo_radius, mars_orbit_radius)
    v_hel_mars = np.sqrt(mu_sun * (2 / d_MS - 1 / transfer_a))
    v_mars_orbit = np.sqrt(mu_sun / d_MS)
    v_inf_mars = v_hel_mars - v_mars_orbit
    
    r_apo_mars, v_mars_apo, v_mars_per, a_mars = compute_capture_orbit(mars_orbit_radius)

    if not inclination_midcourse and not use_aerobraking:
        delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
        delta_v_inclination = compute_mars_inclination_change(v_mars_circ)
        delta_v_drag = compute_atmospheric_drag(mars_orbit_radius, v_mars_circ, spacecraft_mass, Cd, cross_sectional_area, delta_t)
        print("Aerobraking disabled — inclination change at mars circ.")
    elif not inclination_midcourse and use_aerobraking:
        delta_v_2 = np.abs(v_mars_per - v_mars_circ)
        delta_v_inclination = compute_mars_inclination_change(v_mars_apo)
        delta_v_drag = compute_atmospheric_drag(mars_orbit_radius, v_mars_per, spacecraft_mass, Cd, cross_sectional_area, delta_t)
        print("Aerobraking enabled — inclination change at mars apo.")
    elif inclination_midcourse and use_aerobraking:
        delta_v_2 = np.abs(v_mars_per - v_mars_circ)
        delta_v_inclination = 20 / 1000  # 20 m/s
        delta_v_drag = compute_atmospheric_drag(mars_orbit_radius, v_mars_per, spacecraft_mass, Cd, cross_sectional_area, delta_t)
        print("Inclination changed midcourse — aerobraking enabled.")
    else:
        delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
        delta_v_inclination = 20 / 1000
        delta_v_drag = compute_atmospheric_drag(mars_orbit_radius, v_mars_circ, spacecraft_mass, Cd, cross_sectional_area, delta_t)
        print("Inclination changed midcourse — aerobraking disabled.")

    delta_v_station_keeping = compute_station_keeping(years=4.5, delta_v_per_year=60)
    delta_v_deorbit = compute_deorbit()

    total_delta_v_launcher = (
        #delta_v_launch
        delta_v_1+delta_v_inclination
    )

    total_delta_v_spacecraft = (delta_v_2 + delta_v_deorbit)
    period = compute_mars_period(mars_orbit_radius)
    print(f"Period of Mars orbit: {period/60:.2f} minutes")
    # === OUTPUT ===
    print(f"\nEccentricity of transfer orbit: {e:.4f}")
    print("=== Mars Mission ΔV Budget ===")
    #print(f"ΔV: Launch to LEO              = {delta_v_launch:.3f} km/s")
    print(f"ΔV: Mars Transfer Injection      = {delta_v_1:.3f} km/s")
    print(f"ΔV: Inclination Change         = {delta_v_inclination:.3f} km/s")
    print(f"ΔV: Mars Capture & Circularization       = {delta_v_2:.3f} km/s")
    print(f"ΔV: Station Keeping (4.5 yrs)  = {delta_v_station_keeping:.3f} km/s")
    print(f"ΔV: End-of-Life Deorbit        = {delta_v_deorbit:.3f} km/s")
    print(f"ΔV: Atmospheric Drag           = {delta_v_drag:.3f} km/s")
    print("------------------------------------------")
    print(f"Total Mission ΔV (launcher)    = {total_delta_v_launcher:.3f} km/s")
    print(f"Total Mission ΔV (spacecraft)  = {total_delta_v_spacecraft:.3f} km/s")
    return {
        "Mars Transfer Injection": delta_v_1,
        "Inclination Change": delta_v_inclination,
        "Capture & Circularization": delta_v_2,
        "Station Keeping": delta_v_station_keeping,
        "Deorbit": delta_v_deorbit,
        "Atmospheric Drag": delta_v_drag
    }

def run_all_scenarios():
    scenarios = {
        "Aerobraking + Midcourse":      dict(use_aerobraking=True,  inclination_midcourse=True),
        "No Aerobraking + Midcourse":   dict(use_aerobraking=False, inclination_midcourse=True),
        "Aerobraking + No Midcourse":   dict(use_aerobraking=True,  inclination_midcourse=False),
        "No Aerobraking + No Midcourse":dict(use_aerobraking=False, inclination_midcourse=False),
    }
    results = {}
    for label, opts in scenarios.items():
        dv = main(
            use_aerobraking=opts["use_aerobraking"],
            inclination_midcourse=opts["inclination_midcourse"],
            leo_alt=200, mars_orbit_alt=200,
            spacecraft_mass=850, Cd=2.6, cross_sectional_area=2,
            delta_t=3.3*365*24*3600
        )
        results[label] = dv
    return results


#SIMULATION

def plot_comparison_deltav(results):
    mission_phases = [
        "Mars Transfer Injection",
        "Inclination Change",
        "Capture & Circularization",
        "Station Keeping",
        "Deorbit",
        "Atmospheric Drag"
    ]
    plt.figure(figsize=(10, 6))
    for label, dv_dict in results.items():
        dv_list = [dv_dict[phase] for phase in mission_phases]
        cumulative = np.cumsum(dv_list)
        plt.plot(mission_phases, cumulative, marker='o', label=label)
    plt.title("Cumulative ΔV Budget per Mission Phase")
    plt.xlabel("Mission Phase")
    plt.ylabel("Cumulative ΔV (km/s)")
    plt.xticks(rotation=30)
    plt.grid(False)
    plt.legend()
    plt.tight_layout()
    plt.show()


def simulate_and_plot():
    # Define ranges for LEO and Mars orbit altitudes
    leo_altitudes = np.arange(150, 500, 50)  # LEO altitudes in km
    mars_altitudes = np.arange(200, 1000, 100)  # Mars orbit altitudes in km

    # Prepare data storage for plots
    scenarios = {
        "With Aerobraking, No Midcourse": {"use_aerobraking": True, "inclination_midcourse": False, "results": []},
        "Without Aerobraking, No Midcourse": {"use_aerobraking": False, "inclination_midcourse": False, "results": []},
        "With Aerobraking, With Midcourse": {"use_aerobraking": True, "inclination_midcourse": True, "results": []},
        "Without Aerobraking, With Midcourse": {"use_aerobraking": False, "inclination_midcourse": True, "results": []},
    }

    # Iterate over scenarios
    for scenario_name, scenario_params in scenarios.items():
        results = []
        for leo_alt in leo_altitudes:
            for mars_alt in mars_altitudes:
                # Calculate total ΔV for the given scenario and altitudes
                leo_radius = r_earth + leo_alt
                mars_orbit_radius = r_mars + mars_alt
                v_leo = np.sqrt(mu_earth / leo_radius)
                v_mars_circ = np.sqrt(mu_mars / mars_orbit_radius)

                delta_v_launch = compute_launch_to_leo(leo_radius)
                #delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)
                delta_v_1, transfer_a, e, v_inf_leo = compute_dir_tranfer_injection(leo_radius, mars_orbit_radius)

                v_hel_mars = np.sqrt(mu_sun * (2 / d_MS - 1 / transfer_a))
                v_mars_orbit = np.sqrt(mu_sun / d_MS)
                v_inf_mars = v_hel_mars - v_mars_orbit

                r_apo_mars, v_mars_apo, v_mars_per, a_mars = compute_capture_orbit(mars_orbit_radius)

                if not scenario_params["inclination_midcourse"] and not scenario_params["use_aerobraking"]:
                    delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
                    delta_v_inclination = compute_mars_inclination_change(v_mars_circ)
                elif not scenario_params["inclination_midcourse"] and scenario_params["use_aerobraking"]:
                    delta_v_2 = np.abs(v_mars_per - v_mars_circ)
                    delta_v_inclination = compute_mars_inclination_change(v_mars_apo)
                elif scenario_params["inclination_midcourse"]:
                    delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
                    delta_v_inclination = 20 / 1000  # 20 m/s
                else:
                    delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
                    delta_v_inclination = 20 / 1000

                delta_v_station_keeping = compute_station_keeping(years=4.5, delta_v_per_year=60)
                delta_v_deorbit = compute_deorbit()

                total_delta_v = (
                    delta_v_launch +
                    delta_v_1 +
                    delta_v_inclination +
                    delta_v_2 +
                    delta_v_station_keeping +
                    delta_v_deorbit
                )

                # Store results
                results.append((leo_alt, mars_alt, total_delta_v))

        # Store results for the scenario
        scenario_params["results"] = results

    # Plot results
    for scenario_name, scenario_params in scenarios.items():
        results = scenario_params["results"]
        leo_altitudes = [r[0] for r in results]
        mars_altitudes = [r[1] for r in results]
        delta_vs = [r[2] for r in results]

        plt.figure(figsize=(10, 6))
        sc = plt.scatter(leo_altitudes, mars_altitudes, c=delta_vs, cmap='viridis')
        plt.colorbar(sc, label="Total ΔV (km/s)")
        plt.title(f"ΔV Budget: {scenario_name}")
        plt.xlabel("LEO Altitude (km)")
        plt.ylabel("Mars Orbit Altitude (km)")
        plt.grid()
        plt.show()
    

def simulate_and_plot2():
    # Define ranges for LEO and Mars orbit altitudes
    leo_altitudes = np.arange(150, 500, 50)  # LEO altitudes in km
    mars_altitudes = np.arange(200, 1000, 100)  # Mars orbit altitudes in km

    # Prepare data storage for plots
    scenarios = {
        "With Aerobraking, No Midcourse": {"use_aerobraking": True, "inclination_midcourse": False, "results": []},
        "Without Aerobraking, No Midcourse": {"use_aerobraking": False, "inclination_midcourse": False, "results": []},
        "With Aerobraking, With Midcourse": {"use_aerobraking": True, "inclination_midcourse": True, "results": []},
        "Without Aerobraking, With Midcourse": {"use_aerobraking": False, "inclination_midcourse": True, "results": []},
    }

    # Iterate over scenarios
    for scenario_name, scenario_params in scenarios.items():
        results = []
        for leo_alt in leo_altitudes:
            for mars_alt in mars_altitudes:
                # Calculate total ΔV for the given scenario and altitudes
                leo_radius = r_earth + leo_alt
                mars_orbit_radius = r_mars + mars_alt
                v_leo = np.sqrt(mu_earth / leo_radius)
                v_mars_circ = np.sqrt(mu_mars / mars_orbit_radius)

                delta_v_launch = compute_launch_to_leo(leo_radius)
                #delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)
                delta_v_1, transfer_a, e, v_inf_leo = compute_dir_tranfer_injection(leo_radius, mars_orbit_radius)
                v_hel_mars = np.sqrt(mu_sun * (2 / d_MS - 1 / transfer_a))
                v_mars_orbit = np.sqrt(mu_sun / d_MS)
                v_inf_mars = v_hel_mars - v_mars_orbit

                r_apo_mars, v_mars_apo, v_mars_per, a_mars = compute_capture_orbit(mars_orbit_radius)

                if not scenario_params["inclination_midcourse"] and not scenario_params["use_aerobraking"]:
                    delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
                    delta_v_inclination = compute_mars_inclination_change(v_mars_circ)
                elif not scenario_params["inclination_midcourse"] and scenario_params["use_aerobraking"]:
                    delta_v_2 = np.abs(v_mars_per - v_mars_circ)
                    delta_v_inclination = compute_mars_inclination_change(v_mars_apo)
                elif scenario_params["inclination_midcourse"] and not scenario_params["use_aerobraking"]:
                    delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
                    delta_v_inclination = 20 / 1000  # 20 m/s
                else:
                    delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
                    delta_v_inclination = 20 / 1000

                delta_v_station_keeping = compute_station_keeping(years=4.5, delta_v_per_year=60)
                delta_v_deorbit = compute_deorbit()

                total_delta_v = (
                    delta_v_launch +
                    delta_v_1 +
                    delta_v_inclination +
                    delta_v_2 +
                    delta_v_station_keeping +
                    delta_v_deorbit
                )

                # Store results
                results.append((leo_alt, mars_alt, total_delta_v))

        # Store results for the scenario
        scenario_params["results"] = results

    # Plot results
    for scenario_name, scenario_params in scenarios.items():
        results = scenario_params["results"]

        # Extract data for LEO altitude and Mars orbit altitude
        leo_altitudes = sorted(set(r[0] for r in results))
        mars_altitudes = sorted(set(r[1] for r in results))

        # Calculate average ΔV for each LEO altitude and Mars orbit altitude
        delta_v_vs_leo = [np.mean([r[2] for r in results if r[0] == leo_alt]) for leo_alt in leo_altitudes]
        delta_v_vs_mars = [np.mean([r[2] for r in results if r[1] == mars_alt]) for mars_alt in mars_altitudes]

        # Create subplots
        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        # Plot ΔV vs LEO altitude
        axs[0].plot(leo_altitudes, delta_v_vs_leo, marker='o', label="ΔV vs LEO Altitude")
        axs[0].set_title(f"ΔV vs LEO Altitude ({scenario_name})")
        axs[0].set_xlabel("LEO Altitude (km)")
        axs[0].set_ylabel("Total ΔV (km/s)")
        axs[0].grid()
        axs[0].legend()

        # Plot ΔV vs Mars orbit altitude
        axs[1].plot(mars_altitudes, delta_v_vs_mars, marker='o', label="ΔV vs Mars Orbit Altitude", color='orange')
        axs[1].set_title(f"ΔV vs Mars Orbit Altitude ({scenario_name})")
        axs[1].set_xlabel("Mars Orbit Altitude (km)")
        axs[1].set_ylabel("Total ΔV (km/s)")
        axs[1].grid()
        axs[1].legend()

        # Show the plots
        plt.tight_layout()
        plt.show()


# Run the simulation and plot


if __name__ == "__main__":
    main(use_aerobraking=True, inclination_midcourse= True, leo_alt=200, mars_orbit_alt = 120, spacecraft_mass = 1044.6, Cd = 2.705, cross_sectional_area = 1.7*1.2, delta_t = 1.88*365*24*3600)
    # main(use_aerobraking=False, inclination_midcourse= False, leo_alt=200, mars_orbit_alt = 212, spacecraft_mass = 700, Cd = 2.6, cross_sectional_area = 2, delta_t = 3.8*365*24*3600)
    # main(use_aerobraking=True, inclination_midcourse= False, leo_alt=200, mars_orbit_alt = 212, spacecraft_mass = 700, Cd = 2.6, cross_sectional_area = 2, delta_t = 3.3*365*24*3600)
    # main(use_aerobraking=False, inclination_midcourse= True, leo_alt=200, mars_orbit_alt = 212, spacecraft_mass = 700, Cd = 2.6, cross_sectional_area = 2, delta_t = 3.8*365*24*3600)
    period = compute_mars_period(212.48+r_mars)
    print(f"Period of Mars orbit: {period/60:.2f} minutes")
    # simulate_and_plot()
    # simulate_and_plot2()
    #plot_comparison_deltav(results=run_all_scenarios())
