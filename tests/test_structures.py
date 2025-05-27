import pytest
import numpy as np
import matplotlib.pyplot as plt
import Structures.Structures2 as st

def test_Polygon_element_lengths():
    ### INPUTS ###
    x = [0, 1, 1, 0] #x-coordinates of the polygon points
    y = [0, 0, 1, 1] #y-coordinates of the polygon points
    stringers = [0, 0, 0, 0] #number of stringers per element
    t = np.linspace(0.0001, 0.006, 50) #thickness (range) of the skin

    test_box = st.Polygon(height=4.5, px=x, py=y, stringers=stringers, thickness=t)

    expected_lengths = [1, 1, 1, 1]

    assert test_box.elements == expected_lengths


def test_Stringer_hat_stress():
    aluminium = st.Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = st.Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)

    assert np.isclose(hat_stringer.crip_stress, 316.84e6, atol=1e5)


def test_Stringer_hat_area():
    aluminium = st.Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = st.Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)

    expected_area = 2.6e-5

    assert np.isclose(hat_stringer.area, expected_area, atol=1e-6)


def test_Load_calculation_launcher_loading():
    ### INPUTS ###
    x = [0, 2, 2, 0] #x-coordinates of the polygon points
    y = [0, 0, 2, 2] #y-coordinates of the polygon points
    stringers = [2, 2, 2, 2] #number of stringers per element
    t = 0.001 #thickness (range) of the skin
    sc_mass = 691 #mass of the spacecraft in kg
    sc_height = 1 #height of the spacecraft in m

    ### CREATE OBJECTS ###
    aluminium = st.Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = st.Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)
    box = st.Polygon(height=sc_height, px=x, py=y, stringers=stringers, thickness=t)
    load_calc = st.Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass)

    assert np.isclose(load_calc.launcher_loading(), 60461.5, atol=1e-1) 


def test_Load_calculation_calculate_poly_loads():
    ### INPUTS ###
    x = [0, 2, 2, 0] #x-coordinates of the polygon points
    y = [0, 0, 2, 2] #y-coordinates of the polygon points
    stringers = [2, 2, 2, 2] #number of stringers per element
    t = np.array([0.001]) #thickness (range) of the skin
    sc_mass = 691 #mass of the spacecraft in kg
    sc_height = 1 #height of the spacecraft in m

    ### CREATE OBJECTS ###
    aluminium = st.Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = st.Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)
    box = st.Polygon(height=sc_height, px=x, py=y, stringers=stringers, thickness=t)
    load_calc = st.Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass)

    expected_stress = 17446247.95
    expected_area = 8.208e-3
    expected_load = expected_stress * expected_area
    load = load_calc.calculate_poly_loads()
    
    assert np.isclose(load, expected_load, atol=1e4)
    assert np.isclose(load_calc.total_stress, expected_stress, atol=1e4)
    assert np.isclose(load_calc.total_area, expected_area, atol=1e-5)


def test_Load_calculation_calculate_weight():
    ### INPUTS ###
    x = [0, 2, 2, 0] #x-coordinates of the polygon points
    y = [0, 0, 2, 2] #y-coordinates of the polygon points
    stringers = [2, 2, 2, 2] #number of stringers per element
    t = np.array([0.001]) #thickness (range) of the skin
    sc_mass = 691 #mass of the spacecraft in kg
    sc_height = 1 #height of the spacecraft in m

    ### CREATE OBJECTS ###
    aluminium = st.Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = st.Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)
    box = st.Polygon(height=sc_height, px=x, py=y, stringers=stringers, thickness=t)
    load_calc = st.Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass)

    expected_weight = 8.208e-3 * 2800

    load = load_calc.calculate_poly_loads()
    weight = load_calc.calculate_weight()

    assert np.isclose(weight, expected_weight, atol=1e-5)


def test_Load_calculation_calculate_area_inertia_thickness():
    ### INPUTS ###
    x = [0, 2, 2, 0] #x-coordinates of the polygon points
    y = [0, 0, 2, 2] #y-coordinates of the polygon points
    stringers = [2, 2, 2, 2] #number of stringers per element
    t = np.array([0.001]) #thickness (range) of the skin
    sc_mass = 691 #mass of the spacecraft in kg
    sc_height = 1 #height of the spacecraft in m

    ### CREATE OBJECTS ###
    aluminium = st.Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = st.Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)
    box = st.Polygon(height=sc_height, px=x, py=y, stringers=stringers, thickness=t)
    load_calc = st.Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass)

    launch_load = load_calc.launcher_loading()
    rigidity_t = load_calc.calculate_area_inertia_thickness()
    
    expected_rigidity_t = 1.8556e-5

    assert np.isclose(rigidity_t, expected_rigidity_t, atol=1e-7)
