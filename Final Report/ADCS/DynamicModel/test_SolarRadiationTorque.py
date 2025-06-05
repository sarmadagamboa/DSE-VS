from SolarRadiationTorque import maximum_torque
from ModelConstants import *

def test_solar_radiation_torque():
    torque = maximum_torque(q, A_s, c_ps, c_g_s)
    expected_torque = np.array([-0.000098034, 0.000000000000000, 0.000000000000000])
    assert np.allclose(torque, expected_torque, atol=1e-12), f"Expected {expected_torque}, but got {torque}"