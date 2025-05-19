import pytest
import numpy as np
from Astrodynamics import mars_delta_v_budget_functions as m

def test_compute_launch_to_leo():
    leo_radius = m.r_earth + 200  # km
    v = m.compute_launch_to_leo(leo_radius)
    assert isinstance(v, float)
    assert v > 0

def test_compute_transfer_injection():
    leo_radius = m.r_earth + 200
    mars_orbit_radius = m.r_mars + 200
    v_leo = m.compute_launch_to_leo(leo_radius)
    result = m.compute_transfer_injection(leo_radius, mars_orbit_radius, v_leo)
    assert isinstance(result, tuple)
    assert len(result) == 4
    for val in result:
        assert isinstance(val, float)

def test_compute_dir_tranfer_injection():
    r_earth = m.r_earth + 200
    mars_orbit_radius = m.r_mars + 200
    result = m.compute_dir_tranfer_injection(r_earth, mars_orbit_radius)
    assert isinstance(result, tuple)
    assert len(result) == 4
    for val in result:
        assert isinstance(val, float)

def test_compute_capture_orbit():
    mars_orbit_radius = m.r_mars + 200
    result = m.compute_capture_orbit(mars_orbit_radius)
    assert isinstance(result, tuple)
    assert len(result) == 4
    for val in result:
        assert isinstance(val, float)

def test_compute_mars_inclination_change():
    v_orbit = 3.5  # km/s
    dv = m.compute_mars_inclination_change(v_orbit)
    assert isinstance(dv, float)
    assert dv > 0

def test_compute_mars_circularization():
    v_inf_mars = 2.5
    mars_orbit_radius = m.r_mars + 200
    v_mars_circ = 3.5
    dv = m.compute_mars_circularization(v_inf_mars, mars_orbit_radius, v_mars_circ)
    assert isinstance(dv, float)

def test_compute_mars_period():
    mars_orbit_radius = m.r_mars + 200
    T = m.compute_mars_period(mars_orbit_radius)
    assert isinstance(T, float)
    assert T > 0

def test_compute_station_keeping():
    dv = m.compute_station_keeping(years=4.5, delta_v_per_year=60)
    assert isinstance(dv, float)
    assert dv > 0

def test_compute_deorbit():
    dv = m.compute_deorbit()
    assert isinstance(dv, float)
    assert dv > 0

def test_compute_atmospheric_drag():
    dv = m.compute_atmospheric_drag(
        mars_altitude=400,
        spacecraft_velocity=3500,
        spacecraft_mass=1000,
        drag_coefficient=2.2,
        cross_sectional_area=2.0,
        delta_t=1000
    )
    assert isinstance(dv, float)

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

def test_plot_comparison_deltav_runs(monkeypatch):
    # Just check that the plotting function runs without error
    results = m.run_all_scenarios()
    # Patch plt.show to avoid blocking
    monkeypatch.setattr(m.plt, "show", lambda: None)
    m.plot_comparison_deltav(results)