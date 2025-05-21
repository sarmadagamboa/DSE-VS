import numpy as np
import pytest
import Thermal_Control as tc


def test_filter_subsystems_above_Tcold():
    """Subsystems whose Tmin is higher than the threshold should be returned."""
    T_cold = -50.0
    cold_list = [-20.0, 0.0, -60.0, -80.0, 15.0, 5.0]
    names = ["A", "B", "C", "D", "E", "F"]

    expected = [("A", -20.0), ("B", 0.0), ("E", 15.0), ("F", 5.0)]  # hand calc
    assert tc.filter_subsystems_above_Tcold(T_cold, cold_list, names) == expected


def test_compute_payload_heater_power():
    """Payload heater power for T_cold = 10 C."""
    T_cold = 10.0
    expected = 18.759  # hand calc
    assert tc.compute_payload_heater_power(T_cold) == pytest.approx(expected, rel=1e-6)

def test_compute_satellite_area():
    """Surface area of a 2x3x4 m rectangular prism."""
    expected = 52.0  # hand calc
    assert tc.compute_satellite_area(2.0, 3.0, 4.0) == pytest.approx(expected)

def test_compute_TCS_mass():
    """TCS mass with A_rad = 2 m2, A_tot = 10 m2."""
    expected = 6.38  # hand calc
    assert tc.compute_TCS_mass(2.0, 10.0) == pytest.approx(expected)

def test_compute_fluxes():
    """Direct, albedo, and IR fluxes for unity view factor."""
    q_sol, q_alb, q_ir = tc.compute_fluxes(view_factor=1.0)
    assert q_sol == pytest.approx(586.2)          # hand calc
    assert q_alb == pytest.approx(146.55)         # hand calc
    assert q_ir  == pytest.approx(109.851248)     # hand calc

def test_compute_radiator_area():
    """Radiator area when external loads are zero."""
    expected = 1.831552014  # hand calc
    assert tc.compute_radiator_area(0.0, 0.0, 0.0) == pytest.approx(expected, rel=1e-6)

def test_compute_cold_temperature():
    """Equilibrium cold-case temperature for A_rad = 2 m2, no external fluxes."""
    expected = 5.254267  # hand calc (deg C)
    temp = tc.compute_cold_temperature(2.0, 0.0, 0.0, 0.0)
    assert temp == pytest.approx(expected, abs=1e-3)

def test_compute_heater_power_coldcase_zero():
    """A small radiator should need no extra heater power."""
    power = tc.compute_heater_power_coldcase(0.0, 0.0, 0.0, 1.0)
    assert power == 0.0

def test_compute_heater_power_coldcase_positive():
    """A large radiator should require positive heater power."""
    expected = 1005.826146  # hand calc
    power = tc.compute_heater_power_coldcase(0.0, 0.0, 0.0, 5.0)
    assert power == pytest.approx(expected, rel=1e-6)

def test_main(capsys):
    """main() should run without error and print a key headline."""
    tc.main()
    out = capsys.readouterr().out
    assert "Radiator area" in out
