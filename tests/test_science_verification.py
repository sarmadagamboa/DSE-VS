import pytest
import numpy as np
import matplotlib.pyplot as plt
from unittest.mock import patch, MagicMock
import sys
import os


sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

@pytest.fixture
def mock_show():
    with patch('matplotlib.pyplot.show') as mock_show:
        yield mock_show


@pytest.fixture
def mock_plt():
    with patch('matplotlib.pyplot', autospec=True) as mock_plt:
        yield mock_plt


def test_calculation_formulas(): # Test formulas used in each measurement technique
    
    altitude = 200000  # m
    R = 3.3895e6  # Mars radius (m)
    mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
    r = R + altitude
    n = np.sqrt(mu / r**3)
    t_3yr_earth = 3 * 365 * 24 * 60 * 60  # 3 Earth years in seconds
    
    # Test l values
    l_range = np.array([1, 10, 100, 400], dtype=np.float64)
    
    # Test for LRI+ACC
    epsilon = 1e-8
    f = 1
    N = f * t_3yr_earth
    d = 100000
    gamma = 2 * np.arctan((d / 2) / r)
    
    Fl = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
    Fl *= (1 / np.abs(np.sin(l_range * gamma / 2)))
    
    sigma_lri_acc = (epsilon / np.sqrt(N)) * (1 / (n * r)) * ((r / R)**l_range) * Fl
    
    # Test for LRI+CAI
    epsilon = 1e-10
    d = 80000
    gamma = 2 * np.arctan((d / 2) / r)
    
    Fl = l_range / (np.sqrt(2 + 4 * l_range + 4 * l_range**2))
    Fl *= (1 / np.abs(np.sin(l_range * gamma / 2)))
    
    sigma_lri_cai = (epsilon / np.sqrt(N)) * (1 / (n * r)) * ((r / R)**l_range) * Fl
    
    # Test for QGG
    epsilon = 1e-10
    d = 0
    gamma = 2 * np.arctan((d / 2) / r)
    
    Fl = 1 / np.sqrt(2 * (1 + 2*l_range + 2*l_range**2) * (2 + 2*l_range + l_range**2))
    sigma_qgg = (epsilon / np.sqrt(N)) * (1 / n**2) * (r / R)**l_range * Fl
    
    # Test for Doppler
    epsilon = 5e-6
    f = 0.0167
    N = f * t_3yr_earth
    d = 0
    gamma = 2 * np.arctan((d / 2) / r)
    
    Fl = l_range / np.sqrt(1/2 + l_range + l_range**2)
    sigma_doppler = (epsilon / np.sqrt(N)) * (1 / (n * r)) * (r / R)**l_range * Fl
    
    # Test signal power law
    signal_power = 8.5e-5 / l_range**2
    
    # verify calculations
    assert sigma_lri_acc.shape == l_range.shape
    assert sigma_lri_cai.shape == l_range.shape
    assert sigma_qgg.shape == l_range.shape
    assert sigma_doppler.shape == l_range.shape
    
    
    # Verify signal power law
    l_idx = 1  # Index for l=10 in l_range
    expected_signal = 8.5e-5 / 10**2
    assert abs(signal_power[l_idx] - expected_signal) < 1e-15


def test_separations_handling():  # Test separations parameter is handled correctly 
   
    def calculate_with_separation(technique_name, separation):
        
        altitude = 200000  # m
        R = 3.3895e6  # Mars radius (m)
        mu = 4.282837e13  # Mars gravitational parameter (m^3/s^2)
        r = R + altitude
        n = np.sqrt(mu / r**3)
        
        # Test value
        l = 10

        gamma = 2 * np.arctan((separation / 2) / r)
        
        if technique_name in ['LRI+ACC (GRACE-FO)', 'LRI+CAI']:
            Fl = l / (np.sqrt(2 + 4 * l + 4 * l**2))
            Fl *= (1 / np.abs(np.sin(l * gamma / 2)))
            return gamma, Fl
        return gamma, None
    
    separations = {
        'LRI+ACC (GRACE-FO)': 100000,  
        'LRI+CAI': 80000,                              
    }
    
    # Calculate gamma and Fl
    lri_acc_gamma, lri_acc_Fl = calculate_with_separation('LRI+ACC (GRACE-FO)', separations['LRI+ACC (GRACE-FO)'])
    lri_cai_gamma, lri_cai_Fl = calculate_with_separation('LRI+CAI', separations['LRI+CAI'])
    
    # Check  gamma and Fl  values 
    assert lri_acc_gamma > lri_cai_gamma > 0
    assert lri_acc_Fl is not None
    assert lri_cai_Fl is not None