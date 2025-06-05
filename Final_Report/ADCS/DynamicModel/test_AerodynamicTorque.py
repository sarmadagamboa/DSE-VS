from AerodynamicTorque import aerodynamic_torque
from ModelConstants import *

def test_aerodynamic_torque():
    Torque = aerodynamic_torque(C_d, rho, A_aero, r, c_g_a, c_pa)
    expected_torque = np.array([0.00000000, 9.2006011e-6, 0.00000000]) 
    assert np.allclose(Torque, expected_torque, atol=1e-6), f"Expected {expected_torque}, but got {Torque}"