import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from ModelConstants import *
from SolarRadiationTorque import maximum_torque
from AerodynamicTorque import aerodynamic_torque
from GravityTorque import gravity_torque

# Inertia matrix
I = np.array([[Ixx, -Ixy, -Ixz],
              [-Ixy, Iyy, -Iyz],
              [-Ixz, -Iyz, Izz]])
inv_I = np.linalg.inv(I)

def disturbance(t):
    Solar_pressure_torque = maximum_torque(q, A_s, c_ps, c_g_s)
    Gravity_gradient_torque = gravity_torque(r=r, theta=theta, I_y=Iyy, I_z=Izz)
    Aerodynamic_torque = aerodynamic_torque(C_d=C_d, rho=rho, A=A_aero, r=r, c_g=c_g_a, c_pa=c_pa)
    if t < 1000:
        Internal_disturbance_torque = np.array([0, 0, 0])
    elif t >= 1000 and t < 4600:
        Internal_disturbance_torque = np.array([2*1.083e-5, 0, 0]) # moving the antenna
    else:
        Internal_disturbance_torque = np.array([0, 0, 0])
    # if t < 3600:
    #     Internal_disturbance_torque = np.array([0, 0, 0])
    # elif t >= 3600 and t < 3620: 
    #     Internal_disturbance_torque = np.array([1e-4, 1e-4, 1e-4])
    # else:
    #     Internal_disturbance_torque = np.array([np.sin(t/60)*1e-5, np.cos(t/60)*1e-5, np.sin(t/120)*1e-5])
    
    total_torque = Solar_pressure_torque + Gravity_gradient_torque + Aerodynamic_torque + Internal_disturbance_torque
    # total_torque = np.array([0, 0, 0])
    return total_torque

def controller(omega, attitude_error, omega_desired=np.array([om_x, om_y, om_z])):
    Kp = np.array([1, 1, 1])   # Proportional gains for long term stability: [1, 1, 1]
    Kd = np.array([50, 50, 50])      # Derivative gains for long term stability: [10, 15, 10]
    tau = np.zeros(3)  # Initialize torque vector
    tau[0] = -Kp[0] * attitude_error[0] - Kd[0] * (omega[0] - omega_desired[0])
    tau[1] = -Kp[1] * attitude_error[1] - Kd[1] * (omega[1] - omega_desired[1])
    tau[2] = -Kp[2] * attitude_error[2] - Kd[2] * (omega[2] - omega_desired[2])
    # Clamp and quantize torque
    tau = np.clip(tau, min_torque, max_torque)
    tau[0] = min_x_bit * np.round(tau[0] / min_x_bit)
    tau[1] = min_y_bit * np.round(tau[1] / min_y_bit)
    tau[2] = min_z_bit * np.round(tau[2] / min_z_bit)
    return tau

def dynamics(t, state):
    omega = state[0:3]             
    attitude = state[3:6]           

    attitude_error = attitude - desired_attitude(t)  # Error in attitude
    tau = controller(omega, attitude_error)
    d = disturbance(t)

    omega_cross = np.array([
        [0, -omega[2], omega[1]],
        [omega[2], 0, -omega[0]],
        [-omega[1], omega[0], 0]
    ])

    omega_dot = inv_I @ (tau + d - omega_cross @ (I @ omega))
    attitude_dot = omega  # Approximate for small angles

    return np.concatenate([omega_dot, attitude_dot])

def plot(omega_x, omega_y, omega_z, roll, pitch, yaw, time, sol):
    # Desired attitude profiles
    desired_r = np.zeros_like(roll)
    desired_p = -sol.t * 2 * np.pi / (110*60)
    desired_y = np.zeros_like(yaw)

    Settling_time = int(len(time) / 10)  # Assuming the first 10% of the time is settling time
    # Max deviations
    print(f'Maximum deviation from desired roll: {np.max(np.abs(roll[Settling_time:] - desired_r[Settling_time:]))*1000} mrad')
    print(f'Maximum deviation from desired pitch: {np.max(np.abs(pitch[Settling_time:] - desired_p[Settling_time:]))*1000} mrad')
    print(f'Maximum deviation from desired yaw: {np.max(np.abs(yaw[Settling_time:] - desired_y[Settling_time:]))*1000} mrad')

    # Torque logs
    torque_log = []
    for i in range(len(sol.t)):
        omega_i = sol.y[0:3, i]
        attitude_i = sol.y[3:6, i]
        desired_a = desired_attitude(sol.t[i])
        attitude_error = attitude_i - desired_a
        tau_i = controller(omega_i, attitude_error)
        torque_log.append(tau_i)

    torque_log = np.array(torque_log)

    dt = sol.t[1] - sol.t[0]
    print(f"total torque supplied over one orbit: {np.sum(np.abs(torque_log), axis=0) * dt} Nm-s")
    print(f"saturation of moment wheels: {np.sum(torque_log, axis=0)*dt} Nm-s")
    print(f"total disturbance torque over one orbit: {np.sum(np.abs([disturbance(t) for t in sol.t]), axis=0) * dt} Nm-s")
    fig, axs = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

    # --- 1. Angular Velocity ---
    axs[0].plot(sol.t[5:], omega_x[5:], label='Omega_x')
    axs[0].plot(sol.t[5:], omega_y[5:], label='Omega_y')
    axs[0].plot(sol.t[5:], omega_z[5:], label='Omega_z')
    axs[0].set_ylabel('Angular Velocity (rad/s)')
    axs[0].set_title('Spacecraft Angular Velocity')
    axs[0].grid(True)
    axs[0].legend()

    # --- 2. Attitude Angles ---
    axs[1].plot(sol.t[5:], roll[5:], label='Roll')
    axs[1].plot(sol.t[5:], pitch[5:], label='Pitch')
    axs[1].plot(sol.t[5:], yaw[5:], label='Yaw')
    axs[1].set_ylabel('Angle (rad)')
    axs[1].set_title('Spacecraft Attitude Angles')
    axs[1].grid(True)
    axs[1].legend()

    # --- 3. Control Torque ---
    axs[2].plot(sol.t, torque_log[:, 0], label='Torque_x')
    axs[2].plot(sol.t, torque_log[:, 1], label='Torque_y')
    axs[2].plot(sol.t, torque_log[:, 2], label='Torque_z')
    axs[2].set_ylabel('Control Torque (Nm)')
    axs[2].set_title('Control Torque over Time')
    axs[2].grid(True)
    axs[2].legend()

    # --- 4. Disturbance Torque ---
    disturbance_torque = np.array([disturbance(t) for t in sol.t])
    axs[3].plot(sol.t, disturbance_torque[:, 0], label='Disturbance Torque X', color='blue')
    axs[3].plot(sol.t, disturbance_torque[:, 1], label='Disturbance Torque Y', color='orange')
    axs[3].plot(sol.t, disturbance_torque[:, 2], label='Disturbance Torque Z', color='green')
    axs[3].set_ylabel('Disturbance Torque (Nm)')
    axs[3].set_title('Disturbance Torque over Time')
    axs[3].grid(True)
    axs[3].legend()
    # --- Common X-axis ---
    axs[3].set_xlabel('Time (s)')
    # Layout adjustment
    plt.tight_layout()
    plt.show()

def solve(t_span, state0, length):
    t_eval = np.linspace(t_span[0], t_span[1], length)
    sol = solve_ivp(dynamics, t_span, state0, t_eval=t_eval, method='RK45')
    return sol.y[0], sol.y[1], sol.y[2], sol.y[3], sol.y[4], sol.y[5], sol.t, sol

def desired_attitude(t):
    # Desired attitude profile # if i change desired here, I need to also change it in the plot function
    desired_roll = 0
    desired_pitch = -t * 2 * np.pi / (110 * 60)  # One full revolution over one orbit
    desired_yaw = 0
    return np.array([desired_roll, desired_pitch, desired_yaw])


if __name__ == "__main__":
    # Initial state:
    omega0 = np.array([om_x, om_y, om_z])
    attitude0 = np.array([0, 0, 0])
    state0 = np.concatenate([omega0, attitude0])

    # Desired angular velocity
    omega_desired = np.array([om_x, om_y, om_z])

    # Time span and evaluation points
    t_span = (0, 110 * 60) # if orbit length is changed here, it also needs to be changed in desired_attitude and plot
    num_points = 100000  # Number of points for evaluation

    # Solve the ODE
    omega_x, omega_y, omega_z, roll, pitch, yaw, time, sol = solve(t_span, state0, num_points)

    # plotting
    plot(omega_x, omega_y, omega_z, roll, pitch, yaw, time, sol)
    
