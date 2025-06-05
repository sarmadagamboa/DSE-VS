import numpy as np
from GravityTorque import gravity_torque
from ModelConstants import *

def test_gravity_torque():
    Torque = gravity_torque(r, theta, Iyy, Izz)
    expected_torque = np.array([0.00000000, 0.00000000, 0.00000000]) 
    assert np.allclose(Torque, expected_torque, atol=1e-6), f"Expected {expected_torque}, but got {Torque}"
    