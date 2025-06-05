import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from ModelConstants import *
from SolarRadiationTorque import maximum_torque
from AerodynamicTorque import aerodynamic_torque
from GravityTorque import gravity_torque

# inertia matrix
I = np.array([[Ixx, -Ixy, -Ixz],
              [-Ixy, Iyy, -Iyz],
              [-Ixz, -Iyz, Izz]])
inv_I = np.linalg.inv(I)

def disturbance(t):
    Solar_pressure_torque = maximum_torque(q, A_s, c_ps, c_g_s)  # Solar radiation torque
    
    Gravity_gradient_torque = gravity_torque(r=r, theta=theta, I_y=Iyy, I_z=Izz)  # Gravity gradient torque

    Aerodynamic_torque = aerodynamic_torque(C_d=C_d, rho=rho, A=A_s, r=r, c_g=c_g_a, c_pa=c_pa)  # Aerodynamic torque

    if t < 3600 or t >= 3600:  # First hour of the orbit
        Internal_disturbance_torque = np.array([0, 0, 0])  # No internal disturbance torque in the first hour
    # elif t >= 3600 and t < 3620: 
    #     Internal_disturbance_torque = np.array([1, -1, 1])
    # else:
    #     Internal_disturbance_torque = np.array([np.sin(t/60), np.cos(t/60), np.sin(t/120)])

    total_torque = Solar_pressure_torque + Gravity_gradient_torque + Aerodynamic_torque + Internal_disturbance_torque
    return total_torque

def controller(omega):
    Kpx = 100
    Kpy = 100
    Kpz = 100

    tau = omega - omega_desired  # Proportional control to maintain desired angular velocity
    tau[0] = -Kpx * tau[0]  # Control torque in the x direction
    tau[1] = -Kpy * tau[1]  # Control torque in the y direction
    tau[2] = -Kpz * tau[2]  # Control torque in the z direction
    tau = np.clip(tau, min_torque, max_torque)  # Limit the control torque to the specified range

    # Ensure the control torque is at least the minimum torque bits
    tau[0] = min_x_bit * np.round(tau[0] / min_x_bit)
    tau[1] = min_y_bit * np.round(tau[1] / min_y_bit)
    tau[2] = min_z_bit * np.round(tau[2] / min_z_bit)
    return tau


def dynamics(t, omega):
    # Skew-symmetric matrix for omega
    omega_cross = np.array([
        [0, -omega[2], omega[1]],
        [omega[2], 0, -omega[0]],
        [-omega[1], omega[0], 0]
    ])
    
    # Compute control input and disturbance
    tau = controller(omega)
    d = disturbance(t)

    # Compute angular acceleration
    omega_dot = inv_I @ (tau + d - omega_cross @ (I @ omega))
    return omega_dot

# Initial angular velocity
omega0 = np.array([0, 0, 0])  # rad/s

# Desired angular velocity (for control)
omega_desired = np.array([om_x, om_y, om_z])  # rad/s

# Time span and evaluation points
t_span = (0, 110*60)  # one orbbit
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# Integrate
sol = solve_ivp(dynamics, t_span, omega0, t_eval=t_eval, method='RK45')

plt.plot(sol.t[5:], sol.y[0][5:], label='Omega_x')
plt.plot(sol.t[5:], sol.y[1][5:], label='Omega_y')
plt.plot(sol.t[5:], sol.y[2][5:], label='Omega_z')
plt.xlabel('Time (s)')
plt.ylabel('Angular Velocity (rad/s)')
plt.title('Spacecraft Angular Velocity with Feedback Control')
plt.grid(True)
plt.legend()
plt.show()

dt = sol.t[1] - sol.t[0]
omega_x = sol.y[0]
omega_y = sol.y[1]
omega_z = sol.y[2]
roll = []
pitch = []
yaw = []

for i in range(len(omega_x)):
    if i == 0:
        roll.append(0)
        pitch.append(0)
        yaw.append(0)
    else:
        roll.append(roll[i-1]+omega_x[i]*dt)
        pitch.append(pitch[i-1]+omega_y[i]*dt)
        yaw.append(yaw[i-1]+omega_z[i]*dt)
        
plt.plot(sol.t[5:], roll[5:], label='Roll')
plt.plot(sol.t[5:], pitch[5:], label='Pitch')
plt.plot(sol.t[5:], yaw[5:], label='Yaw')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.title('Spacecraft Attitude Angles with Feedback Control')
plt.grid(True)
plt.legend()
plt.show()

desired_roll = np.zeros_like(roll)  # Desired roll angle (constant)
print(f'Maximum deviation from desired roll: {max(abs(roll - desired_roll))} rad')
desired_pitch = -sol.t[5:]*2*np.pi/(110*60)  # Desired pitch angle over one orbit
print(f'Maximum deviation from desired pitch: {max(abs(pitch[5:] - desired_pitch))} rad')
desired_yaw = np.zeros_like(yaw)  # Desired yaw angle (constant)
print(f'Maximum deviation from desired yaw: {max(abs(yaw - desired_yaw))} rad')

torque_log = []
time_log = []

for i in range(len(sol.t)):
    omega_i = sol.y[:, i]
    tau_i = controller(omega_i)
    torque_log.append(tau_i)
    time_log.append(sol.t[i])

torque_log = np.array(torque_log)
time_log = np.array(time_log)
plt.plot(time_log, torque_log[:, 0], label='Torque_x')
plt.plot(time_log, torque_log[:, 1], label='Torque_y')
plt.plot(time_log, torque_log[:, 2], label='Torque_z')
plt.xlabel('Time (s)')
plt.ylabel('Control Torque (Nm)')
plt.title('Control Torque over Time')
plt.grid(True)
plt.legend()
plt.show()

plt.plot(sol.t[5:], abs(pitch[5:] - desired_pitch), label='Pitch Deviation')
plt.xlabel('Time (s)')
plt.ylabel('Deviation (rad)')
plt.title('Pitch Deviation over Time')
plt.grid(True)
plt.legend()
plt.show()