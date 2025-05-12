
import numpy as np
from astropy.constants import G, M_earth, R_earth
from astropy import units as u
import matplotlib.pyplot as plt
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

# Mission parameters


def compute_launch_to_leo(leo_radius):
    return np.sqrt(mu_earth / leo_radius)

def compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo):
    transfer_a = (leo_radius + (d_MS + d_ES) + mars_orbit_radius) / 2
    v_earth_orbit = np.sqrt(mu_sun / d_ES)
    v_hel_LEO = np.sqrt(mu_sun * (2 / d_ES - 1 / transfer_a))
    v_inf_leo = v_hel_LEO - v_earth_orbit
    v_per_LEO = np.sqrt((2 * mu_earth / leo_radius) + v_inf_leo**2)
    delta_v = np.abs(v_per_LEO - v_leo)
    eccentricity = 1 - ((leo_radius + d_ES) / transfer_a)
    return delta_v, transfer_a, eccentricity, v_inf_leo

def compute_capture_orbit(mars_orbit_radius):
    e_mars = 0.8537
    a_mars = mars_orbit_radius / (1 - e_mars)
    r_apo_mars = 2 * a_mars - mars_orbit_radius
    v_mars_apo = np.sqrt(mu_mars / r_apo_mars)
    v_mars_per = np.sqrt(mu_mars * (-1 / a_mars + 2 / mars_orbit_radius))
    return r_apo_mars, v_mars_apo, v_mars_per, a_mars

def compute_mars_inclination_change(v_orbit):

    delta_v_incl = np.sqrt(2 * v_orbit**2 * (1 - np.cos(np.radians(93))))
    return delta_v_incl

def compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ):
    v_apo_mars = np.sqrt((2 * mu_mars / mars_orbit_radius) + v_inf_mars**2)
    return np.abs(v_mars_circ - v_apo_mars)

def compute_station_keeping(years=4.5, delta_v_per_year=60):
    return delta_v_per_year * years / 1000  # km/s

def compute_deorbit():
    return 60 / 1000  # km/s

# # Run all calculations
# e_mars = 0.8537
# a_mars = mars_orbit_radius / (1 - e_mars)
# r_apo_mars = 2 * a_mars - mars_orbit_radius
# v_mars_apo = np.sqrt(mu_mars / r_apo_mars)
# v_mars_per = np.sqrt(mu_mars*(-1/a_mars + 2/mars_orbit_radius))


# delta_v_launch = compute_launch_to_leo()
# delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection()
# v_hel_mars = np.sqrt(mu_sun * (2 / d_MS - 1 / transfer_a))
# v_mars_orbit = np.sqrt(mu_sun / d_MS)
# v_inf_mars = v_hel_mars - v_mars_orbit

# # r_apo_mars = compute_mars_inclination_change(v_inf_mars)

# # delta_v_2 = 0  # Aerobraking assumed
# use_aerobraking = True  
# inclination_midcourse = False
# # if use_aerobraking:
# #     delta_v_2 = 0
# #     print("Aerobraking enabled — circularization ΔV saved.")
# # else:
    
# #     delta_v_inclination = compute_mars_inclination_change(v_mars_circ)[0]
# #     print("Aerobraking disabled — full circularization burn required.")

# if not inclination_midcourse and not use_aerobraking:
#     delta_v_2 = compute_mars_circularization(v_inf_mars)
#     delta_v_inclination = compute_mars_inclination_change(v_mars_circ)[0]
#     print("Aerobraking enabled — circularization ΔV saved.")

# elif not inclination_midcourse and use_aerobraking:
#     delta_v_2 = np.abs(v_mars_per - v_mars_circ)
#     delta_v_inclination = compute_mars_inclination_change(v_mars_apo)[0]
# elif inclination_midcourse and not use_aerobraking:
#     delta_v_2 = compute_mars_circularization(v_inf_mars)
#     delta_v_inclination = 20/1000
# else:
#     delta_v_2 = compute_mars_circularization(v_inf_mars)
#     delta_v_inclination = 20/1000

# delta_v_station_keeping = compute_station_keeping()
# delta_v_deorbit = compute_deorbit()

# # Set to False if no aerobraking

# # Total ΔV
# total_delta_v = (
#     delta_v_launch +
#     delta_v_1 +
#     delta_v_inclination +
#     delta_v_2 +
#     delta_v_station_keeping +
#     delta_v_deorbit
# )

# # Print Results
# print(f'Eccentricity of transfer orbit: {e:.4f}')
# print("=== Mars Mission ΔV Budget ===")
# print(f'ΔV: Launch to LEO              = {delta_v_launch:.3f} km/s')
# print(f'ΔV: LEO to Mars Transfer       = {delta_v_1:.3f} km/s')
# print(f'ΔV: Inclination Change         = {delta_v_inclination:.3f} km/s')
# print(f'ΔV: MARS Circularization       = {delta_v_2:.3f} km/s')
# print(f'ΔV: Station Keeping (4.5 yrs)  = {delta_v_station_keeping:.3f} km/s')
# print(f'ΔV: End-of-Life Deorbit        = {delta_v_deorbit:.3f} km/s')
# print(f'------------------------------------------')
# print(f'Total Mission ΔV               = {total_delta_v:.3f} km/s')


def main(use_aerobraking=True, inclination_midcourse=False, leo_alt=200, mars_orbit_alt=400):
    leo_radius = r_earth + leo_alt
    mars_orbit_radius = r_mars + mars_orbit_alt
    v_leo = np.sqrt(mu_earth / leo_radius)
    v_mars_circ = np.sqrt(mu_mars / mars_orbit_radius)
    
    # === CALCULATIONS ===
    delta_v_launch = compute_launch_to_leo(leo_radius)
    delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)
    
    v_hel_mars = np.sqrt(mu_sun * (2 / d_MS - 1 / transfer_a))
    v_mars_orbit = np.sqrt(mu_sun / d_MS)
    v_inf_mars = v_hel_mars - v_mars_orbit
    
    r_apo_mars, v_mars_apo, v_mars_per, a_mars = compute_capture_orbit(mars_orbit_radius)

    if not inclination_midcourse and not use_aerobraking:
        delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
        delta_v_inclination = compute_mars_inclination_change(v_mars_circ)
        print("Aerobraking disabled — full circularization burn required.")
    elif not inclination_midcourse and use_aerobraking:
        delta_v_2 = np.abs(v_mars_per - v_mars_circ)
        delta_v_inclination = compute_mars_inclination_change(v_mars_apo)
        print("Aerobraking enabled — circularization ΔV saved.")
    elif inclination_midcourse:
        delta_v_2 = compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
        delta_v_inclination = 20 / 1000  # 20 m/s
        print("Inclination changed midcourse.")
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

    # === OUTPUT ===
    print(f"\nEccentricity of transfer orbit: {e:.4f}")
    print("=== Mars Mission ΔV Budget ===")
    print(f"ΔV: Launch to LEO              = {delta_v_launch:.3f} km/s")
    print(f"ΔV: LEO to Mars Transfer       = {delta_v_1:.3f} km/s")
    print(f"ΔV: Inclination Change         = {delta_v_inclination:.3f} km/s")
    print(f"ΔV: Mars Circularization       = {delta_v_2:.3f} km/s")
    print(f"ΔV: Station Keeping (4.5 yrs)  = {delta_v_station_keeping:.3f} km/s")
    print(f"ΔV: End-of-Life Deorbit        = {delta_v_deorbit:.3f} km/s")
    print("------------------------------------------")
    print(f"Total Mission ΔV               = {total_delta_v:.3f} km/s")





#SIMULATION


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
                delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)

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
                delta_v_1, transfer_a, e, v_inf_leo = compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)

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
    main(use_aerobraking=True, inclination_midcourse=False, leo_alt=200, mars_orbit_alt=400)
    simulate_and_plot()
    #simulate_and_plot2()