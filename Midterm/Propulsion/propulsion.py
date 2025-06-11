import numpy as np
import math

#constants
G = 9.80665 #m/s^2
m_dry = 580 #kg

# Mission parameters
mission_duration = 4.5  # years
t_transfer = 8  # months
t_station = 20*290.1 + 365*24*60/10#365*24*60/10#45  # minutes #3.3 years
t_capture = 45 # minutes #30 for biprop or monopropellant and 15*24*60 for biprop/electric

#SST LRI-ACC--> biprop
#SST LRI-CAI and QGG --> biprop + electric
#DT --> monopropellant

#Delta-V requirements
deltav_station = 0.125*290.1++196.21 #0.105#196.21#0.105*70+196.21 #0.285 #60 * (mission_duration - t_transfer/12) #m/s station keeping (60m/s per year): 2nd iteration: 270
deltav_capture = 37#1249  #m/s capture (Assuming aerobreaking assist) --> second iteration: 1000 to 1249

class PropulsionProperties:
    #Option 1a - Electric Propulsion
    ELECTRIC_MAIN = {'name': 'Solar Electric Propulsion',
        'isp': 1700,  # specific impulse (s)
        'eta': 0.5,  # efficiency
        'thruster_mass': 12.5,  # kg (BHT-6000 HIM of 300mN)
        'propellant_density': 1350,  # kg/m^3 (Xenon at room temp, high pressure)
        'feed_system_factor': 0.95,  # 5% mass for feed system
    }
    # Option 1b - Electric Propulsion
    ELECTRIC_SK = {'name': 'Solar Electric Propulsion',
        'isp': 1700,  # specific impulse (s)
        'eta': 0.5,  # efficiency
        'thruster_mass': 2,  # kg (BHT-6000 HIM of 300mN)
        'propellant_density': 1350,  # kg/m^3 (Xenon at room temp, high pressure)
        'feed_system_factor': 0.95,  # 5% mass for feed system
    }
    #Option 2 - Bipropellant Propulsion
    BIPROP = {'name': 'Bipropellant Propulsion',
        'isp': 300,  # specific impulse (s)
        'thruster_mass': 5,  # kg
        'fuel_density': 875,  # kg/m^3 (MMH)
        'oxidizer_density': 1442,  # kg/m^3 (N2O4)
        'pressurant_density': 125,  # kg/m^3 (Helium)
        'ox_fuel_ratio': 0.85,
        'press_to_prop': 7/925,
        'feed_system_factor': 0.80,  # 20% mass for feed system
    }

    BIPROP_SK = {'name': 'Bipropellant Propulsion'}
    #Option 3 - Hybrid Propulsion
    COLDGAS = {
        'name': 'Biprop_SK',
        'isp': 300,  # specific impulse for cold gas (s)
        'thruster_mass': 0.1,  # kg
        'propellant_density': 1200,  # kg/m^3 (N2)
        'feed_system_factor':0.20,  # 20% mass for feed system
    }

    MONOPROP = {
        'name': 'Monopropellant Propulsion',
        'isp': 230,  # specific impulse for cold gas (s)
        'thruster_mass': 4,  # kg
        'propellant_density': 1010,  # kg/m^3 (Hydrazine)
        'feed_system_factor': 0.90,  # 10% mass for feed system
    }

#structural parameters
tank_thickness = 0.003 #m
tank_material_density = 4430  # kg/m^3 (Titanium)

def compute_mass_after_deltav (m_final, deltav, Isp, g=G):
    m_initial = m_final * math.exp(deltav / (Isp * g))
    m_prop = m_initial - m_final
    return m_initial, m_prop

def compute_mass_fuel(m_prop, ox_fuel_ratio,press_to_prop):
    fuel_mass = m_prop/(1+ox_fuel_ratio)
    oxidiser_mass =ox_fuel_ratio*fuel_mass
    pressurant_mass = (fuel_mass+oxidiser_mass)*press_to_prop
    return fuel_mass, oxidiser_mass, pressurant_mass

def compute_tank_volume(m_prop, density):
    return m_prop/density

def compute_spherical_tank_mass(volume, thickness= tank_thickness, rho= tank_material_density):
    tank_radius = ((3 * volume) / (4 * math.pi)) ** (1 / 3)
    surface_area = 4 * math.pi * tank_radius ** 2
    mass = surface_area * thickness * rho
    return tank_radius, mass

def compute_exhaust_velocity (Isp, g = G):
    return Isp*g

def compute_average_thrust(m_initial, deltav, duration_minutes):
    duration_s = duration_minutes * 60
    thrust_avg = m_initial*deltav/duration_s
    return thrust_avg

def compute_power_required(thrust,ve,efficiency):
    return thrust*ve/(2*efficiency)

def analyse_electric_option():
    prop = PropulsionProperties.ELECTRIC_MAIN

    m_before_stkeeping, m_prop_stkeeping = compute_mass_after_deltav(m_dry, deltav_station, prop['isp'])
    m_before_capture, m_prop_capture = compute_mass_after_deltav(m_before_stkeeping, deltav_capture, prop['isp'])

    V_stkeeping = compute_tank_volume(m_prop_stkeeping, prop['propellant_density'])
    V_capture = compute_tank_volume(m_prop_capture, prop['propellant_density'])
    V_total = V_stkeeping + V_capture


    tank_radius, tank_mass = compute_spherical_tank_mass(V_total)

    ve = compute_exhaust_velocity(prop['isp'])

    thrust_stkeeping = compute_average_thrust(m_before_stkeeping, deltav_station, t_station)
    thrust_capture = compute_average_thrust(m_before_capture, deltav_capture, t_capture)

    power_stkeeping = compute_power_required(thrust_stkeeping, ve, prop['eta'])
    power_capture = compute_power_required(thrust_capture, ve, prop['eta'])

    total_mass_propulsion = (2 * prop['thruster_mass'] + tank_mass) / prop['feed_system_factor']

    return {'name': prop['name'],
        'launch_mass': m_before_capture,
        'propellant_mass': {
            'station_keeping': m_prop_stkeeping,
            'capture': m_prop_capture,
            'total': m_prop_stkeeping + m_prop_capture
            },
        'tank_volumes': {
            'total': V_total
        },
        'tank_radii': {
            'propellant': tank_radius
        },
        'thrust': {
            'station_keeping': thrust_stkeeping,
            'capture': thrust_capture
        },
        'power_required': {
            'station_keeping': power_stkeeping,
            'capture': power_capture,
            'total': power_stkeeping + power_capture
        },
        'propulsion_system_mass': total_mass_propulsion
    }

def analyse_biprop_option():
    prop = PropulsionProperties.BIPROP

    m_before_stkeeping, m_prop_stkeeping = compute_mass_after_deltav(m_dry, deltav_station, prop['isp'])
    m_before_capture, m_prop_capture = compute_mass_after_deltav(m_before_stkeeping, deltav_capture, prop['isp'])

    fuel_mass, oxidizer_mass, pressurant_mass = compute_mass_fuel(
        m_prop_capture + m_prop_stkeeping,
        prop['ox_fuel_ratio'],
        prop['press_to_prop']
    )

    total_propellant_mass = fuel_mass + oxidizer_mass + pressurant_mass

    v_fuel = compute_tank_volume(fuel_mass, prop['fuel_density'])
    v_oxidizer = compute_tank_volume(oxidizer_mass, prop['oxidizer_density'])
    v_pressurant = compute_tank_volume(pressurant_mass, prop['pressurant_density'])

    tank_radius_fuel, tank_mass_fuel = compute_spherical_tank_mass(v_fuel)
    tank_radius_oxidizer, tank_mass_oxidizer = compute_spherical_tank_mass(v_oxidizer)
    tank_radius_pressurant, tank_mass_pressurant = compute_spherical_tank_mass(v_pressurant)

    ve = compute_exhaust_velocity(prop['isp'])

    thrust_stkeeping = compute_average_thrust(m_before_stkeeping, deltav_station, t_station)
    thrust_capture = compute_average_thrust(m_before_capture, deltav_capture, t_capture)

    total_mass_propulsion = (prop['thruster_mass'] + tank_mass_fuel +
                             tank_mass_pressurant + tank_mass_oxidizer) / prop['feed_system_factor']

    return {'name': prop['name'],
        'launch_mass': m_before_capture,
        'propellant_mass': {
            'station_keeping': m_prop_stkeeping,
            'capture': m_prop_capture,
            'fuel': fuel_mass,
            'oxidizer': oxidizer_mass,
            'pressurant': pressurant_mass,
            'total': total_propellant_mass
        },
        'tank_volumes':{
            'fuel': v_fuel,
            'oxidizer': v_oxidizer,
            'pressurant': v_pressurant
        },
        'tank_radii': {
            'fuel': tank_radius_fuel,
            'oxidizer': tank_radius_oxidizer,
            'pressurant': tank_radius_pressurant
        },
        'thrust': {
            'station_keeping': thrust_stkeeping,
            'capture': thrust_capture
        },
        'propulsion_system_mass': total_mass_propulsion
    }

def analyse_hybrid_option_biprop_coldgas():
    biprop = PropulsionProperties.BIPROP
    coldgas = PropulsionProperties.BIPROP_SK

    # Calculate masses for station keeping using cold gas
    m_before_stkeeping, m_prop_stkeeping_coldgas = compute_mass_after_deltav(m_dry, deltav_station, biprop['isp'])

    # Calculate masses for capture using biprop
    m_before_capture, m_prop_capture = compute_mass_after_deltav(m_before_stkeeping, deltav_capture, biprop['isp'])

    # Calculate biprop component masses (for capture only)
    fuel_mass, oxidizer_mass, pressurant_mass = compute_mass_fuel(
        m_prop_capture,
        biprop['ox_fuel_ratio'],
        biprop['press_to_prop']
    )

    # Calculate tank volumes
    v_coldgas = compute_tank_volume(m_prop_stkeeping_coldgas, coldgas['propellant_density'])
    v_fuel = compute_tank_volume(fuel_mass, biprop['fuel_density'])
    v_oxidizer = compute_tank_volume(oxidizer_mass, biprop['oxidizer_density'])
    v_pressurant = compute_tank_volume(pressurant_mass, biprop['pressurant_density'])

    # Calculate tank sizes and masses
    tank_radius_coldgas, tank_mass_coldgas = compute_spherical_tank_mass(v_coldgas)
    tank_radius_fuel, tank_mass_fuel = compute_spherical_tank_mass(v_fuel)
    tank_radius_oxidizer, tank_mass_oxidizer = compute_spherical_tank_mass(v_oxidizer)
    tank_radius_pressurant, tank_mass_pressurant = compute_spherical_tank_mass(v_pressurant)

    # Calculate thrust
    ve_coldgas = compute_exhaust_velocity(coldgas['isp'])
    ve_biprop = compute_exhaust_velocity(biprop['isp'])
    thrust_stkeeping = compute_average_thrust(m_before_stkeeping, deltav_station, t_station)
    thrust_capture = compute_average_thrust(m_before_capture, deltav_capture, t_capture)

    # Calculate total propulsion system mass (assuming 6 cold gas thrusters for attitude control)
    coldgas_mass = (tank_mass_coldgas + 6 * coldgas['thruster_mass']) / coldgas['feed_system_factor']
    biprop_mass = (2*biprop['thruster_mass'] + tank_mass_fuel + tank_mass_oxidizer + tank_mass_pressurant) / biprop[
        'feed_system_factor']

    total_mass_propulsion = coldgas_mass + biprop_mass

    return {
        'name': 'Bipropellant / Cold Gas Hybrid',
        'launch_mass': m_before_capture,
        'propellant_mass': {
            'cold_gas': m_prop_stkeeping_coldgas,
            'biprop': m_prop_capture,
            'total': m_prop_stkeeping_coldgas + fuel_mass + oxidizer_mass + pressurant_mass
        },
        'tank_volumes':{
            'cold_gas': v_coldgas,
            'fuel': v_fuel,
            'oxidizer': v_oxidizer,
            'pressurant': v_pressurant
        },
        'tank_radii': {
            'cold_gas': tank_radius_coldgas,
            'fuel': tank_radius_fuel,
            'oxidizer': tank_radius_oxidizer,
            'pressurant': tank_radius_pressurant
        },
        'thrust': {
            'station_keeping': thrust_stkeeping,
            'capture': thrust_capture
        },
        'propulsion_system_mass': total_mass_propulsion
    }

def analyse_hybrid_option_biprop_electric():
    biprop = PropulsionProperties.BIPROP
    electric = PropulsionProperties.ELECTRIC_SK

    # Calculate masses for station keeping using cold gas
    m_before_stkeeping, m_prop_stkeeping_electric = compute_mass_after_deltav(m_dry, deltav_station, electric['isp'])

    # Calculate masses for capture using biprop
    m_before_capture, m_prop_capture = compute_mass_after_deltav(m_before_stkeeping, deltav_capture, biprop['isp'])

    # Calculate biprop component masses (for capture only)
    fuel_mass, oxidizer_mass, pressurant_mass = compute_mass_fuel(
        m_prop_capture,
        biprop['ox_fuel_ratio'],
        biprop['press_to_prop']
    )

    # Calculate tank volumes
    v_electric = compute_tank_volume(m_prop_stkeeping_electric, electric['propellant_density'])
    v_fuel = compute_tank_volume(fuel_mass, biprop['fuel_density'])
    v_oxidizer = compute_tank_volume(oxidizer_mass, biprop['oxidizer_density'])
    v_pressurant = compute_tank_volume(pressurant_mass, biprop['pressurant_density'])

    # Calculate tank sizes and masses
    tank_radius_electric, tank_mass_electric = compute_spherical_tank_mass(v_electric)
    tank_radius_fuel, tank_mass_fuel = compute_spherical_tank_mass(v_fuel)
    tank_radius_oxidizer, tank_mass_oxidizer = compute_spherical_tank_mass(v_oxidizer)
    tank_radius_pressurant, tank_mass_pressurant = compute_spherical_tank_mass(v_pressurant)

    # Calculate thrust
    ve_electric = compute_exhaust_velocity(electric['isp'])
    ve_biprop = compute_exhaust_velocity(biprop['isp'])
    thrust_stkeeping = compute_average_thrust(m_before_stkeeping, deltav_station, t_station)
    thrust_capture = compute_average_thrust(m_before_capture, deltav_capture, t_capture)

    power_stkeeping = compute_power_required(thrust_stkeeping, ve_electric, electric['eta'])

    # Calculate total propulsion system mass (assuming 6 cold gas thrusters for attitude control)
    electric_mass = (tank_mass_electric + 2 * electric['thruster_mass']) / electric['feed_system_factor']
    biprop_mass = (2*biprop['thruster_mass'] + tank_mass_fuel + tank_mass_oxidizer + tank_mass_pressurant) / biprop[
        'feed_system_factor']

    total_mass_propulsion = electric_mass + biprop_mass

    return {
        'name': 'Bipropellant / Electric Hybrid',
        'propellant_mass': {
            'electric': m_prop_stkeeping_electric,
            'biprop': m_prop_capture,
            'total': m_prop_stkeeping_electric + m_prop_capture
        },
        'tank_volumes':{
            'electric': v_electric,
            'fuel': v_fuel,
            'oxidizer': v_oxidizer,
            'pressurant': v_pressurant
        },
        'tank_radii': {
            'electric': tank_radius_electric,
            'fuel': tank_radius_fuel,
            'oxidizer': tank_radius_oxidizer,
            'pressurant': tank_radius_pressurant
        },
        'thrust': {
            'station_keeping': thrust_stkeeping,
            'capture': thrust_capture
        },
        'power_required': {
            'station_keeping': power_stkeeping
        },
        'propulsion_system_mass': total_mass_propulsion,
        'electric_system_mass': electric_mass,
        'biprop_system_mass': biprop_mass
    }

def analyse_monoprop_option():
    prop = PropulsionProperties.MONOPROP

    m_before_stkeeping, m_prop_stkeeping = compute_mass_after_deltav(m_dry, deltav_station, prop['isp'])
    m_before_capture, m_prop_capture = compute_mass_after_deltav(m_before_stkeeping, deltav_capture, prop['isp'])

    V_stkeeping = compute_tank_volume(m_prop_stkeeping, prop['propellant_density'])
    V_capture = compute_tank_volume(m_prop_capture, prop['propellant_density'])
    V_total = V_stkeeping + V_capture


    tank_radius, tank_mass = compute_spherical_tank_mass(V_total)

    ve = compute_exhaust_velocity(prop['isp'])

    thrust_stkeeping = compute_average_thrust(m_before_stkeeping, deltav_station, t_station)
    thrust_capture = compute_average_thrust(m_before_capture, deltav_capture, t_capture)

    total_mass_propulsion = (2 * prop['thruster_mass'] + tank_mass) / prop['feed_system_factor']

    return {'name': prop['name'],
        'launch_mass': m_before_capture,
        'propellant_mass': {
            'station_keeping': m_prop_stkeeping,
            'capture': m_prop_capture,
            'total': m_prop_stkeeping + m_prop_capture
            },
        'tank_volumes': {
            'total': V_total
        },
        'tank_radii': {
            'propellant': tank_radius
        },
        'thrust': {
            'station_keeping': thrust_stkeeping,
            'capture': thrust_capture
        },
        'propulsion_system_mass': total_mass_propulsion
    }

if __name__ == "__main__":
    print(analyse_hybrid_option_biprop_electric())

