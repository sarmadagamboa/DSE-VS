from dataclasses import dataclass
import numpy as np
import pytest
from Power import Mission, kelly_cos, solar_array_sizing, battery_sizing

# Kelly–cosine approximation cases
@dataclass
class KellyCase:
    theta: float
    expected: float

KELLY_CASES = [
    KellyCase(0,   1.000),
    KellyCase(30,  0.866),
    KellyCase(50,  0.635),
    KellyCase(60,  0.450),
    KellyCase(80,  0.100),
    KellyCase(85,  0.000),
]

@pytest.mark.parametrize("case", KELLY_CASES)
def test_kelly_cos_tabulated(case: KellyCase):
    """kelly_cos reproduces tabulated values to <0.1 %."""
    assert np.isclose(kelly_cos(case.theta), case.expected, rtol=1e-3)

def test_kelly_cos_linear_interpolation():
    """
    kelly_cos linearly interpolates between the tabulated values.
    Between 30 deg (0.866) and 50 deg (0.635), at 40 deg it should be exactly halfway.
    """
    low = KELLY_CASES[1].expected   # at 30 deg
    high = KELLY_CASES[2].expected  # at 50 deg
    expected = low + (40 - 30) / (50 - 30) * (high - low)
    val_40 = kelly_cos(40)
    assert np.isclose(val_40, expected, rtol=1e-8)

# Solar‐array sizing
def test_solar_array_sizing_baseline():
    """
    For the body‐fixed array at theta=30 deg,
    verify key outputs.
    """
    m = Mission()
    area, p_eol, p_bol, p_sa, inc = solar_array_sizing(m)

    # EOL power
    assert np.isclose(p_eol, 836.4,        rtol=1e-3)
    # BOL power density
    assert np.isclose(p_bol, 129.450546,   rtol=1e-3)
    # Required power in sunlight
    assert np.isclose(p_sa, 1841.846154,   rtol=1e-3)
    # Array area
    assert np.isclose(area, 14.886448,     rtol=1e-3)
    # Incidence for body‐fixed at 30 degrees
    assert np.isclose(inc,  0.866,         rtol=1e-3)

# battery sizing
def test_battery_sizing_baseline():
    """
    Baseline battery sizing using eclipse=2611.68 s.
    """
    m = Mission()
    bat_mass, bat_vol = battery_sizing(m)

    # Battery mass (kg)
    assert np.isclose(bat_mass, 2.8, rtol=1e-3)
    # Battery volume (L)
    assert np.isclose(bat_vol,  2.1, rtol=1e-3)

def test_battery_sizing_eclipse_scaling():
    """
    Doubling eclipse duration should ~double
    the required battery mass,
    """
    m_short = Mission()
    m_long  = Mission()
    m_long.eclipse = 2 * m_short.eclipse

    mass_short, _ = battery_sizing(m_short)
    mass_long,  _ = battery_sizing(m_long)

    assert np.isclose(mass_long / mass_short, 2.0, rtol=2e-2)

def test_integration_power_mass():
    """
    Call the two sizing functions and compute power_mass
    Verifies that all components work together and give the expected power_mass.
    """
    m = Mission()
    area, p_eol, p_bol, p_sa, inc = solar_array_sizing(m)
    bat_mass, bat_volume = battery_sizing(m)
    power_mass = (area * 2.06) + bat_mass + (0.071 * p_eol) + 0.15

    assert np.isclose(power_mass, 93.0, rtol=1e-3)