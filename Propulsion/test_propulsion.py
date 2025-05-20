import numpy as np
import math
import pytest
import propulsion

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