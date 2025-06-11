import numpy as np
import math

#constants
G = 9.80665 #m/s^2
m_dry_with_initial_prop = 639.3 #kg

# Mission parameters
mission_duration = 4.5  # years
t_transfer = 8  # months
t_station = 45  #3.3 years
t_capture = 45 # minutes #30 for biprop

#Delta-V requirements
deltav_station = 121.3 + 60 #orbit maintainance/drag and EOL maneuver
deltav_capture = 1249  #m/s capture (Assuming aerobreaking assist)

class PropulsionProperties:
    #Electric Propulsion - Stationkeeping
    ELECTRIC_SK = {'name': 'Solar Electric Propulsion',
        'isp': 1242.9,  # specific impulse (s) at 188km: 1242.9  at 177km:1409.6
        'power': 73.4,  #system power (W) at 188km:73.4 at 177km: 113.9
        'thruster_mass': 2.3,  # kg ITA
        'thrust': 1.3,  #mN
        'PXFA_mass' : 4,  #2-5kg,
        'PPU_mass':17.5,  #kg
        'propellant_density': 1830,  # kg/m^3 (Xenon at room temp, high pressure) https://electricrocket.org/IEPC/IEPC1991-107.pdf
        'prop_margin': 0.20  # %
                   }

    #Bipropellant Propulsion - Capture
    BIPROP = {'name': 'Bipropellant Propulsion',
        'thrust': 425, #thrust (N)   340-450N
        'D_throat': 0.01645, #m^2
        'exp_area': 330, #nozzle expansion ratio (Area)
        'chamber_pr': 1035000, #chamber pressure (Pa)
        'isp': 321,  # specific impulse (s)
        'fuel_density': 880,  # kg/m^3 (MMH)
        'oxidizer_density': 1370,  # kg/m^3 (MON-3)
        'press_density': 6.87,  # kg/m^3 (Helium)
        'blow_down_ratio':3.5, #commmon for blow-down systems
        'ox_fuel_ratio': 1.65,
        'prop_margin': 0.10, # %
        'final_press':2000000,#Pa
        'storage_temp':300,#K
        'storage_press':20000000, #Pa
        'gas_const': 2077, # J/kg*K

         # Component masses (kg)
        'thruster_mass': 4.3,  # kg
        'prop_fill_drain_valve': 0.09, #kg,
        'press_fill_drain_valve': 0.06, #kg,
        'parallel_check_valve': 0.04, #kg
        'press_regulator': 5.9, #kg
        'pyro_valve': 0.16,#kg
        'iso_valve': 0.545, #kg
        'filter': 0.567, #kg,
        'press_sens': 0.25, #kg

        # Quantities
        'n_thrusters': 2,
        'n_prop_fill_drain_valve': 4,
        'n_press_fill_drain_valve': 11,
        'n_parallel_check_valve': 4,
        'n_press_regulator': 4,
        'n_pyro_valve': 12,
        'n_iso_valve': 18,
        'n_filter': 5,
        'n_press_sens': 8
              }


#structural parameters
tank_thickness = 0.003 #m
tank_material_density = 4430 # kg/m^3 (Titanium)

def thrust_and_characteristic_velocity(Thrust, p_chamb, D_throat, Isp, g):
    A_throat = np.pi*(D_throat/2)**2
    Cf = Thrust / (p_chamb * A_throat)
    c_theor = Isp * g / Cf
    return Cf, c_theor

def compute_mass_after_deltav (m_final, deltav, Isp, g=G):
    m_initial = m_final * math.exp(deltav / (Isp * g))
    m_prop = m_initial - m_final
    return m_initial, m_prop

def compute_mass_fuel(m_prop, ox_fuel_ratio):
    fuel_mass = m_prop/(1+ox_fuel_ratio)
    oxidiser_mass =ox_fuel_ratio*fuel_mass
    return fuel_mass, oxidiser_mass

def compute_tank_volume(m_prop, density):
    return m_prop/density

def compute_pressurant(final_press,storage_press, storage_temp, gas_const, v_prop_biprop):
    M_press =final_press * v_prop_biprop/(gas_const*storage_temp)
    v_press = M_press*gas_const*storage_temp/storage_press
    return M_press, v_press
"""
def compute_pressurant(volume, B,density):
    V_ullage = volume/(B-1)
    ullage_volume_m3 = V_ullage/ 1000  # Convert to m³
    M_press = V_ullage * density
    return V_ullage, M_press
"""
def compute_spherical_tank_mass(volume, thickness= tank_thickness, rho= tank_material_density):
    tank_radius = ((3 * volume) / (4 * math.pi)) ** (1 / 3)
    surface_area = 4 * math.pi * tank_radius ** 2
    volume_L = volume*1000
    mass = 0.0334*volume_L + 10.408
    return tank_radius, mass
def compute_exhaust_velocity (Isp, g = G):
    return Isp*g

biprop = PropulsionProperties.BIPROP
electric = PropulsionProperties.ELECTRIC_SK

thrust_coefficient, characteristic_velocity=thrust_and_characteristic_velocity(biprop['thrust'], biprop['chamber_pr'], biprop['D_throat'], biprop['isp'], g=G)

# Initial propulsion system masses (excluding propellants, calculated later)
# Electric: 2 thrusters + 2* feedsystem + 2*PPU (1 for redundancy)
electric_base_mass = 2 * electric['thruster_mass'] + 2* electric['PXFA_mass'] +2*electric['PPU_mass']
# Biprop: 3 thrusters (1 for redundancy) + 15 fill/drain + 12 Press sens + 5 filters + 4 reg + 5 check + 14 pyro +20 iso
biprop_base_mass = (
    biprop['n_thrusters'] * biprop['thruster_mass'] +
    biprop['n_prop_fill_drain_valve'] * biprop['prop_fill_drain_valve'] +
    biprop['n_press_fill_drain_valve'] * biprop['press_fill_drain_valve'] +
    biprop['n_press_sens'] * biprop['press_sens'] +
    biprop['n_parallel_check_valve'] * biprop['parallel_check_valve'] +
    biprop['n_pyro_valve'] * biprop['pyro_valve'] +
    biprop['n_iso_valve'] * biprop['iso_valve'] +
    biprop['n_filter'] * biprop['filter'] +
    biprop['n_press_regulator'] * biprop['press_regulator']
)

print(f"Electric propulsion base mass: {electric_base_mass:.1f} kg")
print(f"Biprop propulsion base mass: {biprop_base_mass:.1f} kg")

actual_propulsion_hardware = electric_base_mass + biprop_base_mass
print(f"Total actual propulsion hardware(no tanks): {actual_propulsion_hardware:.1f} kg")

#initial_prop_estimate = 70  # kg
#m_dry_corrected = m_dry_with_initial_prop - initial_prop_estimate + actual_propulsion_hardware
m_dry_corrected = 693.3
print(f"\nMass correction:")
#print(f"Original dry mass (with {initial_prop_estimate}kg prop estimate): {m_dry_with_initial_prop:.1f} kg")
print(f"Corrected dry mass (with detailed prop): {m_dry_corrected:.1f} kg")
#print(f"Propulsion mass change: {actual_propulsion_hardware - initial_prop_estimate:+.1f} kg")

#Mass before capture maneuver
m_before_stkeeping, m_prop_stkeeping_electric_nomargin = compute_mass_after_deltav(m_dry_corrected, deltav_station, electric['isp'])
m_prop_stkeeping_electric = m_prop_stkeeping_electric_nomargin + m_prop_stkeeping_electric_nomargin* electric['prop_margin']
print(f"\nStation keeping calculation:")
print(f"Electric propellant needed (no margin): {m_prop_stkeeping_electric_nomargin:.1f} kg")
print(f"Electric propellant needed: {m_prop_stkeeping_electric:.1f} kg")
print(f"Mass before capture (after station keeping): {m_before_stkeeping:.1f} kg")

# Calculate masses for capture using biprop
m_initial_mission, m_prop_capture_nomargin = compute_mass_after_deltav(m_before_stkeeping, deltav_capture, biprop['isp'])
m_prop_capture = m_prop_capture_nomargin + m_prop_capture_nomargin * biprop['prop_margin']
print(f"\nCapture maneuver calculation:")
print(f"Biprop propellant needed: {m_prop_capture:.1f} kg")
print(f"Total initial mission mass: {m_initial_mission:.1f} kg")

# Calculate biprop component masses (for capture only)
fuel_mass, oxidizer_mass = compute_mass_fuel(m_prop_capture,biprop['ox_fuel_ratio'])

# Calculate tank volumes
v_electric = compute_tank_volume(m_prop_stkeeping_electric, electric['propellant_density'])
v_fuel = compute_tank_volume(fuel_mass, biprop['fuel_density'])
v_oxidizer = compute_tank_volume(oxidizer_mass, biprop['oxidizer_density'])
v_prop_biprop = v_fuel + v_oxidizer

# Calculate pressurant requirements
press_mass, v_pressurant = compute_pressurant(biprop['final_press'],biprop ['storage_press'], biprop['storage_temp'], biprop['gas_const'],v_prop_biprop)


print(f"\nBiprop breakdown:")
print(f"\nBiprop breakdown:")
print(f"Fuel mass: {fuel_mass:.1f} kg")
print(f"Oxidizer mass: {oxidizer_mass:.1f} kg")
print(f"Pressurant mass: {press_mass:.1f} kg")

print(f"\nTank volumes:")
print(f"Electric propellant tank: {v_electric*1000:.1f} L")
print(f"Fuel tank: {v_fuel*1000:.1f} L")
print(f"Oxidizer tank: {v_oxidizer*1000:.1f} L")
print(f"Pressurant tank: {v_pressurant*1000:.1f} L")

# Calculate tank sizes and masses
tank_radius_electric, tank_mass_electric = compute_spherical_tank_mass(v_electric)
tank_radius_fuel, tank_mass_fuel = compute_spherical_tank_mass(v_fuel)
tank_radius_oxidizer, tank_mass_oxidizer = compute_spherical_tank_mass(v_oxidizer)
tank_radius_pressurant, tank_mass_pressurant = compute_spherical_tank_mass(v_pressurant)

print(f"\nTank masses:")
print(f"Electric tank mass: {tank_mass_electric:.1f} kg")
print(f"Fuel tank mass: {tank_mass_fuel:.1f} kg")
print(f"Oxidizer tank mass: {tank_mass_oxidizer:.1f} kg")
print(f"Pressurant tank mass: {tank_mass_pressurant:.1f} kg")

# Calculate total system masses (including propellants)
electric_total_mass = tank_mass_electric + electric_base_mass + m_prop_stkeeping_electric
biprop_total_mass = tank_mass_fuel + tank_mass_oxidizer + tank_mass_pressurant +biprop_base_mass + m_prop_capture + press_mass

total_propulsion_mass = electric_total_mass + biprop_total_mass

print(f"\nTotal system masses:")
print(f"Electric system total: {electric_total_mass:.1f} kg")
print(f"Biprop system total: {biprop_total_mass:.1f} kg")
print(f"Total propulsion mass: {total_propulsion_mass:.1f} kg")

 #The issue is we need to account for tank masses in the dry mass

# Tank masses that should be included in dry mass
total_tank_mass = tank_mass_electric + tank_mass_fuel + tank_mass_oxidizer + tank_mass_pressurant

# The dry mass used in rocket equation should include tanks
m_dry_with_tanks = m_dry_corrected + total_tank_mass

propulsion_system_mass = total_tank_mass + actual_propulsion_hardware

print(f"\nCORRECTED CALCULATION:")
print(f"Dry mass with propulsion hardware (no tanks): {m_dry_corrected:.1f} kg")
print(f"Tank structures mass: {total_tank_mass:.1f} kg")
print(f"Total propulsion system mass: {propulsion_system_mass:.1f}kg")
print(f"Dry mass with propulsion and tanks: {m_dry_with_tanks:.1f} kg")

# Recalculate with correct dry mass including tanks
m_before_stkeeping_corrected, m_prop_stkeeping_corrected = compute_mass_after_deltav(
    m_dry_with_tanks, deltav_station, electric['isp']
)

m_initial_corrected, m_prop_capture_corrected = compute_mass_after_deltav(
    m_before_stkeeping_corrected, deltav_capture, biprop['isp']
)

print(f"\nRecalculated with tank-inclusive dry mass:")
print(f"Electric propellant needed: {m_prop_stkeeping_corrected:.1f} kg")
print(f"Biprop propellant needed: {m_prop_capture_corrected:.1f} kg")
print(f"Corrected total mission mass: {m_initial_corrected:.1f} kg")

# Now verify
propellants_only = m_prop_stkeeping_corrected + m_prop_capture_corrected + press_mass
calculated_total_final = m_dry_with_tanks + propellants_only

print(f"\nFINAL VERIFICATION:")
print(f"Dry mass (including tanks): {m_dry_with_tanks:.1f} kg")
print(f"Propellants only: {propellants_only:.1f} kg")
print(f"Final calculated total: {calculated_total_final:.1f} kg")
print(f"Should match mission mass: {m_initial_corrected:.1f} kg")
print(f"Final difference: {abs(calculated_total_final - m_initial_corrected):.3f} kg")

