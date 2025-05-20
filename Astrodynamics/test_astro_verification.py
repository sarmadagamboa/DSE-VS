import pytest
import numpy as np
import matplotlib.pyplot as plt
from Astrodynamics import mars_delta_v_budget_functions as m
from sso_repeat import repeat_sso, generate_repeat_curves, generate_sso_curve, plot_repeat_curves

def test_compute_launch_to_leo():
    leo_radius = m.r_earth + 200  # km
    v = m.compute_launch_to_leo(leo_radius)
    # Manual expected value (approximate, based on GM/r and typical LEO)
    expected = 7.784  # km/s (typical for 200 km LEO)
    assert v == pytest.approx(expected, rel=1e-2)

def test_compute_deorbit():
    dv = m.compute_deorbit()
    expected = 60  # If your function returns 60 (check your implementation)
    assert dv == pytest.approx(expected, rel=1e-6)

def test_compute_dir_tranfer_injection():
    r_earth = m.r_earth 
    mars_orbit_radius = m.r_mars + 200
    result = m.compute_dir_tranfer_injection(r_earth, mars_orbit_radius)
   
    expected_delta_v = 11.4  # Replace with your manual calculation
    assert result[0] == pytest.approx(expected_delta_v, rel=1e-1)

def test_compute_capture_orbit():
    mars_orbit_radius = m.r_mars + 200
    result = m.compute_capture_orbit(mars_orbit_radius)

    expected_r_apo_mars = 45480.9  # km
    expected_v_mars_apo = 0.371  # km/s (approximate)
    expected_v_mars_per = 4.7  # km/s (approximate)
    expected_a_mars = 24535.2  # km (approximate)
    assert result[0] == pytest.approx(expected_r_apo_mars, rel=1e-1)
    assert result[1] == pytest.approx(expected_v_mars_apo, rel=1e-1)
    assert result[2] == pytest.approx(expected_v_mars_per, rel=1e-1)
    assert result[3] == pytest.approx(expected_a_mars, rel=1e-1)

def test_compute_mars_inclination_change():
    inclination = 93  # degrees

    result = m.compute_mars_inclination_change(v_orbit=7)
    expected_delta_v= 5.01 # manual calculation
    assert expected_delta_v == pytest.approx(result, rel=1e-1)

def test_compute_mars_circularization():
    mars_orbit_radius = m.r_mars + 200
    v_inf_mars = 3.5  # km/s
    v_mars_circ = 3.5
    result = m.compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
    expected_delta_v = 2.509  # Replace with your manual calculation
    assert expected_delta_v == pytest.approx(result, rel=1e-1)

def test_compute_mars_period():
    mars_orbit_radius = m.r_mars + 200
    T = m.compute_mars_period(mars_orbit_radius)
    # Manual expected value (approximate, for 200 km Mars orbit)
    expected = 6529  
    assert T == pytest.approx(expected, rel=1e-1)  

def test_compute_station_keeping():
    dv = m.compute_station_keeping(years=4.5, delta_v_per_year=60)
    expected = 270/1000  # 60*4.5
    assert dv == pytest.approx(expected, rel=1e-6)

def test_compute_atmospheric_drag():

    dv_drag = m.compute_atmospheric_drag(
        mars_altitude=200,
        spacecraft_velocity=5,
        spacecraft_mass=800,
        drag_coefficient=2,
        cross_sectional_area=5,
        delta_t=4*365*24*3600
    )
    expected_delta_v = 0.2
    assert dv_drag == pytest.approx(expected_delta_v, rel=1e-1)

def test_main_and_run_all_scenarios():
    result = m.main(
        use_aerobraking=True,
        inclination_midcourse=True,
        leo_alt=200,
        mars_orbit_alt=200,
        spacecraft_mass=850,
        Cd=2.6,
        cross_sectional_area=2,
        delta_t=3.3*365*24*3600
    )
    assert isinstance(result, dict)
    expected_keys = {
        "Mars Transfer Injection",
        "Inclination Change",
        "Capture & Circularization",
        "Station Keeping",
        "Deorbit",
        "Atmospheric Drag"
    }
    assert set(result.keys()) == expected_keys

    all_results = m.run_all_scenarios()
    assert isinstance(all_results, dict)
    for scenario, dv_dict in all_results.items():
        assert set(dv_dict.keys()) == expected_keys

def test_simulate_and_plot_runs(monkeypatch):
    # Patch plt.show to avoid blocking
    monkeypatch.setattr(plt, "show", lambda: None)
    m.simulate_and_plot()  # Should run without error

def test_simulate_and_plot2_runs(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    m.simulate_and_plot2()  # Should run without error

def test_plot_comparison_deltav_runs(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    results = m.run_all_scenarios()
    m.plot_comparison_deltav(results)  # Should run without error

def test_repeat_sso():

    unique_orbits_df, _ = repeat_sso(sol_range=(1, 1), tol=1e-6, max_iter=100, period_bounds_hr=(1.5, 3.0))

    assert np.isclose(unique_orbits_df["Inclination_deg"], 92.647, rtol=1e-2)
    assert np.isclose(unique_orbits_df["Altitude_km"], 296, rtol=1e-1)

def test_plot_repeat_curves():

    plot_repeat_curves()

    assert True