import numpy as np
import math

#constants
G = 9.80665 #m/s^2

# Mission parameters
mission_duration = 4.5  # years
t_transfer = 8  # months
t_station = 20*290.1 + 365*24*60/10  #3.3 years
t_capture = 45 # minutes #30 for biprop

#Delta-V requirements
deltav_station = 56+196.21  #0.105*70+196.21 #orbit maintainance/drag and EOL maneuver
deltav_capture = 1308.14  #m/s capture (Assuming aerobreaking assist)

class PropulsionProperties:
    #Electric Propulsion - Stationkeeping
    ELECTRIC_SK = {'name': 'Solar Electric Propulsion',
        'isp': 1500,  # specific impulse (s) at 188km: 1242.9  at 177km:1409.6
        'power': 1200, #system power (W) at 188km:73.4 at 177km: 113.9
        'thrust': 1.3, #mN #1-20mN
        'propellant_density': 1350, # kg/m^3 (Xenon at room temp, high pressure) https://electricrocket.org/IEPC/IEPC1991-107.pdf
        'prop_margin':0.20,#%

        # Component masses (kg)
        'PXFA_mass' : 1.1, #2-5kg,
        'PPU_mass':20.5, #kg
        'thruster_mass': 3,  # kg ITA
        'tank_mass': 3.5, #kg
        'tank_volume': 0.006 #m^3 (5L)
        }

    #Bipropellant Propulsion - Capture
    BIPROP = {'name': 'Bipropellant Propulsion',
              'thrust': 425,  #thrust (N)   340-450N
              'D_throat': 0.01645,  #m^2
              'exp_area': 330,  #nozzle expansion ratio (Area)
              'chamber_pr': 1035000,  #chamber pressure (Pa)
              'isp': 321,  # specific impulse (s)
              'fuel_density': 880,  # kg/m^3 (MMH)
              'oxidizer_density': 1370,  # kg/m^3 (MON-3)
              'press_density': 6.87,  # kg/m^3 (Helium)
              'blow_down_ratio':3.5,  #commmon for blow-down systems
              'ox_fuel_ratio': 1.65,
              'prop_margin': 0.15,

              'final_press': 2000000,  # Pa
              'storage_temp': 300,  # K
              'storage_press': 20000000,  # Pa
              'gas_const': 2077,  # J/kg*K

              # Component masses (kg)
              'thruster_mass': 4.3,  # kg
              'prop_fill_drain_valve': 0.09,  #kg,
              'press_fill_drain_valve': 0.06,  #kg,
              'parallel_check_valve': 0.04,  #kg
              'press_regulator': 5.9,  #kg
              'pyro_valve': 0.16,  #kg
              'iso_valve': 0.545,  #kg
              'filter': 0.567,  #kg,
              'press_sens': 0.25,  #kg


              'pressurant_tank_mass': 8,  #kg
              'prop_tank_mass':21, #kg

              'pressurant_tank_volume': 0.040,  # m^3
              'tank_volume':0.198, #m^3

              # Quantities
              'n_thrusters': 2,
              'n_prop_fill_drain_valve': 4,
              'n_press_fill_drain_valve': 9,
              'n_parallel_check_valve': 4,
              'n_press_regulator': 4,
              'n_pyro_valve': 8,
              'n_iso_valve': 6,
              'n_filter': 3,
              'n_press_sens': 8
              }

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


biprop = PropulsionProperties.BIPROP
electric = PropulsionProperties.ELECTRIC_SK

# Initial propulsion system masses (excluding propellants, calculated later)
# Electric: 2 thrusters + 2* feedsystem + 2*PPU  + 2*Xenon Tanks
electric_base_mass = 2 * electric['thruster_mass'] + 2* electric['PXFA_mass'] +2*electric['PPU_mass']+2*electric['tank_mass']
biprop_base_mass = (
    biprop['n_thrusters'] * biprop['thruster_mass'] +
    biprop['n_prop_fill_drain_valve'] * biprop['prop_fill_drain_valve'] +
    biprop['n_press_fill_drain_valve'] * biprop['press_fill_drain_valve'] +
    biprop['n_press_sens'] * biprop['press_sens'] +
    biprop['n_parallel_check_valve'] * biprop['parallel_check_valve'] +
    biprop['n_pyro_valve'] * biprop['pyro_valve'] +
    biprop['n_iso_valve'] * biprop['iso_valve'] +
    biprop['n_filter'] * biprop['filter'] +
    biprop['n_press_regulator'] * biprop['press_regulator'] +
    2 * biprop['prop_tank_mass'] + biprop['pressurant_tank_mass']
)

propulsion_hardware = electric_base_mass + biprop_base_mass
print(f"Total propulsion hardware: {propulsion_hardware:.1f} kg")


m_dry_actual = 677


# Calculate station keeping propellant requirements

m_before_stkeeping, m_prop_stkeeping_electric_nomargin = compute_mass_after_deltav(m_dry_actual, deltav_station, electric['isp'])

m_prop_stkeeping_electric = m_prop_stkeeping_electric_nomargin * (1 + electric['prop_margin'])

print(f"\nStation keeping calculation:")
print(f"Electric propellant needed (no margin): {m_prop_stkeeping_electric_nomargin:.1f} kg")
print(f"Electric propellant needed (with margin): {m_prop_stkeeping_electric:.1f} kg")

# Calculate capture maneuver requirements
m_initial_mission, m_prop_capture_nomargin = compute_mass_after_deltav(
    m_before_stkeeping, deltav_capture, biprop['isp']
)
m_prop_capture = m_prop_capture_nomargin * (1 + biprop['prop_margin'])

print(f"\nCapture maneuver calculation:")
print(f"Biprop propellant needed (no margin): {m_prop_capture_nomargin:.1f} kg")
print(f"Biprop propellant needed (with margin): {m_prop_capture:.1f} kg")

# Calculate biprop component masses
fuel_mass, oxidizer_mass = compute_mass_fuel(m_prop_capture, biprop['ox_fuel_ratio'])

# Calculate tank volumes required vs available
v_fuel_required = compute_tank_volume(fuel_mass, biprop['fuel_density'])
v_oxidizer_required = compute_tank_volume(oxidizer_mass, biprop['oxidizer_density'])

print(f"\nBiprop breakdown:")
print(f"Fuel mass: {fuel_mass:.1f} kg")
print(f"Oxidizer mass: {oxidizer_mass:.1f} kg")

print(f"\nTank volume requirements:")
print(f"Fuel volume required: {v_fuel_required*1000:.1f} L")
print(f"Oxidizer volume required: {v_oxidizer_required*1000:.1f} L")

# Calculate pressurant requirements
v_prop_biprop = v_fuel_required + v_oxidizer_required
press_mass, v_pressurant = compute_pressurant(
    biprop['final_press'], biprop['storage_press'],
    biprop['storage_temp'], biprop['gas_const'], v_prop_biprop
)

print(f"Pressurant mass: {press_mass:.1f} kg")
print(f"Pressurant volume required: {v_pressurant*1000:.1f} L")

# Tank capacity checks for all three biprop tanks
print(f"\n TANK CAPACITY CHECKS")

# Check if electric propellant fits in available tank volume
max_electric_prop = 2 * electric['tank_volume'] * electric['propellant_density']
print(f"Maximum electric propellant capacity: {max_electric_prop:.1f} kg")
if m_prop_stkeeping_electric > max_electric_prop:
    print(f"WARNING: Required electric propellant ({m_prop_stkeeping_electric:.1f} kg) exceeds tank capacity!")
else:
    print(f"Electric propellant fits in tanks with {max_electric_prop - m_prop_stkeeping_electric:.1f} kg margin")


# 1. Fuel tank
max_fuel_capacity = biprop['tank_volume']  # Single tank design
print(f"Maximum fuel tank capacity: {max_fuel_capacity*1000:.1f} L")
if v_fuel_required > max_fuel_capacity:
    print(f"WARNING: Required fuel volume ({v_fuel_required*1000:.1f} L) exceeds tank capacity!")
else:
    print(f"Fuel fits in tank with {(max_fuel_capacity - v_fuel_required)*1000:.1f} L margin")

# 2. Oxidizer check
max_oxidizer_capacity = biprop['tank_volume']  # Single tank design
print(f"Maximum oxidizer tank capacity: {max_oxidizer_capacity*1000:.1f} L")
if v_oxidizer_required > max_oxidizer_capacity:
    print(f"WARNING: Required oxidizer volume ({v_oxidizer_required*1000:.1f} L) exceeds tank capacity!")
else:
    print(f"Oxidizer fits in tank with {(max_oxidizer_capacity - v_oxidizer_required)*1000:.1f} L margin")

# 3. Pressurant check
max_pressurant_capacity = biprop['pressurant_tank_volume']  # Single tank design
print(f"Maximum pressurant tank capacity: {max_pressurant_capacity*1000:.1f} L")
if v_pressurant > max_pressurant_capacity:
    print(f"WARNING: Required pressurant volume ({v_pressurant*1000:.1f} L) exceeds tank capacity!")
else:
    print(f"Pressurant fits in tank with {(max_pressurant_capacity - v_pressurant)*1000:.1f} L margin")



# Calculate total system masses
electric_total_mass = electric_base_mass + m_prop_stkeeping_electric
biprop_total_mass = biprop_base_mass + m_prop_capture + press_mass
total_propulsion_mass = electric_total_mass + biprop_total_mass
total_propellant = m_prop_stkeeping_electric + m_prop_capture + press_mass

print(f"\nTotal system masses:")
print(f"Electric system total: {electric_total_mass:.1f} kg")
print(f"Biprop system total: {biprop_total_mass:.1f} kg")
print(f"Total propulsion mass: {total_propulsion_mass:.1f} kg")
print(f"Total propellant mass: {total_propellant:.1f} kg")

# Final verification

calculated_total_final = m_dry_actual + total_propulsion_mass

