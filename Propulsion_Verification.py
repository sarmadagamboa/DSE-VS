import numpy as np
import math
import pytest

#constants
G = 9.80665 #m/s^2
m_dry = 638 #kg

# Mission parameters
mission_duration = 4.5  # years
t_transfer = 8  # months
t_station = (mission_duration * 12 - t_transfer) * 30 * 24 * 60  # minutes
t_capture = 15*24*60 # minutes #30 for biprop or monopropellant and 15*24*60 for biprop/electric

#SST LRI-ACC--> biprop
#SST LRI-CAI and QGG --> biprop + electric
#DT --> monopropellant

#Delta-V requirements
deltav_station = 270 #60 * (mission_duration - t_transfer/12) #m/s station keeping (60m/s per year): 2nd iteration: 270
deltav_capture = 1249 #m/s capture (Assuming aerobreaking assist) --> second iteration: 1000 to 1249

class PropulsionProperties:
    #Option 1a - Electric Propulsion
    ELECTRIC_MAIN = {'name': 'Solar Electric Propulsion',
        'isp': 2000,  # specific impulse (s)
        'eta': 0.7,  # efficiency
        'thruster_mass': 12.5,  # kg (BHT-6000 HIM of 300mN)
        'propellant_density': 1350,  # kg/m^3 (Xenon at room temp, high pressure)
        'feed_system_factor': 0.95,  # 5% mass for feed system
    }
    # Option 1b - Electric Propulsion
    ELECTRIC_SK = {'name': 'Solar Electric Propulsion',
        'isp': 2000,  # specific impulse (s)
        'eta': 0.7,  # efficiency
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
    #Option 3 - Hybrid Propulsion
    COLDGAS = {
        'name': 'Cold Gas',
        'isp': 60,  # specific impulse for cold gas (s)
        'thruster_mass': 0.1,  # kg
        'propellant_density': 175,  # kg/m^3 (N2)
        'feed_system_factor': 0.80,  # 20% mass for feed system
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
    coldgas = PropulsionProperties.COLDGAS

    # Calculate masses for station keeping using cold gas
    m_before_stkeeping, m_prop_stkeeping_coldgas = compute_mass_after_deltav(m_dry, deltav_station, coldgas['isp'])

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
            'total': m_prop_stkeeping_coldgas + m_prop_capture + fuel_mass + oxidizer_mass + pressurant_mass
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
        'launch_mass': m_before_capture,
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
        'propulsion_system_mass': total_mass_propulsion
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
def test_compute_mass_after_deltav():
# Test 1: Zero delta-v should result in no propellant needed
    m_final = 1000  # kg
    deltav = 0  # m/s
    Isp = 300  # s
    m_initial, m_prop = compute_mass_after_deltav(m_final, deltav, Isp)
    assert m_initial == pytest.approx(m_final)
    assert m_prop == pytest.approx(0)

# Test 2: Typical values
    m_final = 1000  # kg
    deltav = 500  # m/s
    Isp = 300  # s
    m_initial, m_prop = compute_mass_after_deltav(m_final, deltav, Isp)
    assert m_initial > m_final
    # Hand calculation verification
    expected_initial = 1000 * math.exp(500 / (300 * 9.80665))
    expected_prop = expected_initial - 1000
    assert m_initial == pytest.approx(expected_initial)
    assert m_prop == pytest.approx(expected_prop)

# Test 3: High Isp should require less propellant
    m_initial_low_isp, m_prop_low_isp = compute_mass_after_deltav(m_final, deltav, Isp=300)
    m_initial_high_isp, m_prop_high_isp = compute_mass_after_deltav(m_final, deltav, Isp=3000)
    assert m_prop_high_isp < m_prop_low_isp


def test_compute_mass_fuel():
    # Test: Simple values
    m_prop = 100  # kg
    ox_fuel_ratio = 1.0  # 1:1 ratio
    press_to_prop = 0.01  # 1% pressurant

    fuel_mass, oxidiser_mass, pressurant_mass = compute_mass_fuel(m_prop, ox_fuel_ratio, press_to_prop)

    assert fuel_mass == pytest.approx(50)
    assert oxidiser_mass == pytest.approx(50)
    assert pressurant_mass == pytest.approx(1)
    assert fuel_mass + oxidiser_mass + pressurant_mass > m_prop  # Due to pressurant

def test_compute_tank_volume():
    # Test 1: Simple values
    m_prop = 100  # kg
    density = 1000  # kg/m^3
    volume = compute_tank_volume(m_prop, density)
    assert volume == pytest.approx(0.1)  # 100 kg / 1000 kg/m^3 = 0.1 m^3

    # Test 2: Extreme density
    m_prop = 100  # kg
    density = 0  # kg/m^3 (very low)
    with pytest.raises(ZeroDivisionError):
        compute_tank_volume(m_prop, density)

def test_compute_spherical_tank_mass():
    # Test 1: Simple values
    volume = 1  # m^3
    thickness = 0.003  # m
    rho = 4430  # kg/m^3

    radius, mass = compute_spherical_tank_mass(volume, thickness, rho)

    # Manual calculation for verification
    expected_radius = ((3 * volume) / (4 * math.pi)) ** (1 / 3)
    expected_surface_area = 4 * math.pi * expected_radius ** 2
    expected_mass = expected_surface_area * thickness * rho

    assert radius == pytest.approx(expected_radius)
    assert mass == pytest.approx(expected_mass)

    # Test 2: Linear scaling with density and thickness
    _, mass1 = compute_spherical_tank_mass(volume = 1, thickness = 0.001, rho =1000)
    _, mass2 = compute_spherical_tank_mass(volume = 1, thickness = 0.001, rho =2000)
    assert mass2 == pytest.approx(2 * mass1)

    _, mass1 = compute_spherical_tank_mass(volume = 1, thickness = 0.001, rho =1000)
    _, mass2 = compute_spherical_tank_mass(volume = 1, thickness = 0.002, rho =1000)
    assert mass2 == pytest.approx(2 * mass1)

def test_compute_exhaust_velocity():
    # Test 1: Standard values
    Isp = 300  # s
    g = 9.80665  # m/s^2
    ve = compute_exhaust_velocity(Isp, g)
    assert ve == pytest.approx(Isp * g)

def test_compute_average_thrust():
    # Test case 1: Standard values
    m_initial = 1000  # kg
    deltav = 500  # m/s
    duration_minutes = 60  # minutes

    thrust = compute_average_thrust(m_initial, deltav, duration_minutes)

    # Test 1: Manual calculation for verification
    expected_thrust = m_initial * deltav / (duration_minutes * 60)
    assert thrust == pytest.approx(expected_thrust)

    # Test 2: Double duration should result in half of the thrust
    thrust1 = compute_average_thrust(m_initial, deltav, duration_minutes= 60)
    thrust2 = compute_average_thrust(m_initial, deltav, duration_minutes = 120)
    assert thrust2 == pytest.approx(thrust1 / 2)

def test_compute_power_required():
    # Test 1: Standard values
    thrust = 1  # N
    ve = 20000  # m/s (typical for electric propulsion)
    efficiency = 0.7

    power = compute_power_required(thrust, ve, efficiency)
    expected_power = thrust * ve / (2 * efficiency)
    assert power == pytest.approx(expected_power)

    # Test 2: Double efficiency should result in half of the power
    power1 = compute_power_required(thrust, ve, efficiency= 0.5)
    power2 = compute_power_required(thrust, ve, efficiency= 1.0)
    assert power2 == pytest.approx(power1 / 2)

def test_analyse_electric_option_integration():
    result = analyse_electric_option()

    # Test overall properties that should be true
    assert result['propellant_mass']['total'] > 0
    assert result['propellant_mass']['capture'] > result['propellant_mass']['station_keeping']

    # Verify relationships that should hold
    assert result['propellant_mass']['total'] == result['propellant_mass']['station_keeping'] + \
            result['propellant_mass']['capture']
    assert result['power_required']['total'] == result['power_required']['station_keeping'] + \
            result['power_required']['capture']

    # Verify total mass relationships
    assert result['launch_mass'] > m_dry

    # Subsystem mass check
    assert result['propulsion_system_mass'] > 0
    assert result['propulsion_system_mass'] < result['launch_mass']
def test_analyse_biprop_option_integration():
    result = analyse_biprop_option()

    assert result['propellant_mass']['total'] > 0
    assert result['propellant_mass']['capture'] > result['propellant_mass']['station_keeping']

    # Verify relationships that should hold
    assert pytest.approx(result['propellant_mass']['total']) == (
            result['propellant_mass']['fuel'] +
            result['propellant_mass']['oxidizer'] +
            result['propellant_mass']['pressurant'])

    # Subsystem mass check
    assert result['propulsion_system_mass'] > 0
    assert result['propulsion_system_mass'] < result['launch_mass']


def test_analyse_monoprop_option_integration():
    result = analyse_monoprop_option()

    # Test overall properties that should be true
    assert result['propellant_mass']['total'] > 0
    assert result['propellant_mass']['capture'] > result['propellant_mass']['station_keeping']

    # Verify relationships that should hold
    assert result['propellant_mass']['total'] == (
            result['propellant_mass']['station_keeping'] +
            result['propellant_mass']['capture'])

    # Verify total mass relationships
    assert result['launch_mass'] > m_dry
    assert result['propulsion_system_mass'] < result['launch_mass']


def test_analyse_hybrid_option_biprop_coldgas_integration():
    result = analyse_hybrid_option_biprop_coldgas()

    # Test overall properties that should be true
    assert result['propellant_mass']['total'] > 0
    assert result['propellant_mass']['biprop'] > result['propellant_mass']['cold_gas']

    # Verify mass relationships
    assert result['launch_mass'] > m_dry
    assert result['propulsion_system_mass'] < result['launch_mass']


def test_analyse_hybrid_option_biprop_electric_integration():
    result = analyse_hybrid_option_biprop_electric()

    # Test overall properties that should be true
    assert result['propellant_mass']['total'] > 0
    assert result['propellant_mass']['biprop'] > result['propellant_mass']['electric']

    # Verify relationships that should hold
    assert result['propellant_mass']['total'] == (
            result['propellant_mass']['electric'] +
            result['propellant_mass']['biprop'])

    # Verify total mass relationships
    assert result['launch_mass'] > m_dry
    assert result['propulsion_system_mass'] < result['launch_mass']

    # Verify power requirements exist and are positive for electric portion
    assert 'power_required' in result
    assert result['power_required']['station_keeping'] > 0



