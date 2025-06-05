import numpy as np
from DynamicModelV3 import disturbance, controller, desired_attitude
from GravityTorque import gravity_torque
from SolarRadiationTorque import maximum_torque
from AerodynamicTorque import aerodynamic_torque
from ModelConstants import *

def test_disturbance():
    expected = gravity_torque(r=r, theta=theta, I_y=Iyy, I_z=Izz) + \
                maximum_torque(q, A_s, c_ps, c_g_s) + \
                aerodynamic_torque(C_d=C_d, rho=rho, A=A_aero, r=r, c_g=c_g_a, c_pa=c_pa)
    actual = disturbance(0)
    assert np.allclose(actual, expected, atol=1e-6), f"Expected {expected}, but got {actual}"
    expected = gravity_torque(r=r, theta=theta, I_y=Iyy, I_z=Izz) + \
                maximum_torque(q, A_s, c_ps, c_g_s) + \
                aerodynamic_torque(C_d=C_d, rho=rho, A=A_aero, r=r, c_g=c_g_a, c_pa=c_pa) + \
                np.array([2*1.083e-5, 0, 0])  # Internal disturbance torque when moving the antenna
    actual = disturbance(2000)
    assert np.allclose(actual, expected, atol=1e-6), f"Expected {expected}, but got {actual}"

def test_controller():
    omega = np.array([np.pi, np.pi, np.pi])  # Initial angular velocity
    attitude_error = np.array([0, 0, 0])  # No error in attitude
    expected_torque = np.array([-0.08, -0.08, -0.08])  # Expected torque after clamping and quantization
    actual_torque = controller(omega, attitude_error)
    assert np.allclose(actual_torque, expected_torque, atol=1e-6), f"Expected maximum torque of {expected_torque} after clamping, but got {actual_torque}"
    omega = np.array([0, -2*np.pi/(110*60), 0])  # Initial angular velocity is the desired state
    attitude_error = np.array([0, 0, 0])  # No error in attitude
    expected_torque = np.array([0, 0, 0])  # Expected torque is zero
    actual_torque = controller(omega, attitude_error)
    assert np.allclose(actual_torque, expected_torque, atol=1e-6), f"Expected zero torque of {expected_torque}, but got {actual_torque}"
    omega = np.array([0, -2*np.pi/(110*60) - 1e-4, 0])  # Initial angular velocity is close to the desired state
    attitude_error = np.array([0, 0, 0])  # No error in attitude
    expected_torque = np.array([0, 0.005, 0])  # Expected torque
    actual_torque = controller(omega, attitude_error)
    assert np.allclose(actual_torque, expected_torque, atol=1e-6), f"Expected torque of {expected_torque}, but got {actual_torque}"

def test_dynamics():
    """
    Will implement this test later
    """
    pass

def test_desired_attitude():
    t = 0
    expected_attitude = np.array([0, 0, 0])  # Initial attitude
    actual_attitude = desired_attitude(t)
    assert np.allclose(actual_attitude, expected_attitude, atol=1e-6), f"Expected {expected_attitude}, but got {actual_attitude}"