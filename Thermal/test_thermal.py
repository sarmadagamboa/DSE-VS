import numpy as np
import pytest
import Thermal_Control as tc

def test_filter_subsystems_above_Tcold():
    tc.filter_subsystems_above_Tcold

def test_compute_payload_heater_power():
    tc.compute_payload_heater_power

def test_compute_satellite_area():
    tc.compute_satellite_area

def test_compute_TCS_mass():
    tc.compute_TCS_mass

def test_compute_fluxes():
    tc.compute_fluxes

def test_compute_radiator_area():
    tc.compute_radiator_area

def test_compute_cold_temperature():
    tc.compute_cold_temperature

def test_compute_heater_power_coldcase():
    tc.compute_heater_power_coldcase

def test_main():
    tc.main

    return True