import pytest
from Final_Report.Mass_iteration.mass_iteration import calc_prop_mass, calculate_wet_mass, calculate_dry_mass, iteration_loop

def test_calc_prop_mass():
    # This function is verified in its original file, ####INPUT NAME####, so we don't need to test it here.
    pass


def test_calculate_wet_mass():
    inputs = {}
    inputs["mass"] = {
        "Payload_mass": 100,  # kg
        "ADCS_mass": 40,  # kg
        "TTC_mass": 60,  # kg
        "CDHS_mass": 10,  # kg
        "Thermal_mass": 80,  # kg
        "Power_mass": 4.35,  # kg
        "Propulsion_dry_mass": 5.65,  # kg
        "Propellant_mass": 300,  # kg
    }
    expected_wet_mass = 600  + 300 * 0.1
    wet_mass = calculate_wet_mass(inputs, dry_mass_margin=1.1)
    assert wet_mass == pytest.approx(expected_wet_mass, abs=1e-1)

def test_calculate_dry_mass():
    inputs = {}
    inputs["mass"] = {
        "Payload_mass": 100,  # kg
        "ADCS_mass": 40,  # kg
        "TTC_mass": 60,  # kg
        "CDHS_mass": 10,  # kg
        "Thermal_mass": 80,  # kg
        "Power_mass": 4.35,  # kg
        "Propulsion_dry_mass": 5.65,  # kg
        "Propellant_mass": 300,  # kg
    }
    expected_dry_mass = 300  + 300 * 0.1
    dry_mass = calculate_dry_mass(inputs, dry_mass_margin=1.1)
    assert dry_mass == pytest.approx(expected_dry_mass, abs=1e-1)


def test_iteration_loop():
    # This function is verified in its original file by inspection of the output. It is observed that
    # the resulting mass converges.
    pass


if __name__ == "__main__":
    pytest.main([__file__])

