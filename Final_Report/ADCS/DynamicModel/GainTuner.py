import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from ModelConstants import *
from SolarRadiationTorque import maximum_torque
from AerodynamicTorque import aerodynamic_torque
from GravityTorque import gravity_torque
import pandas as pd

# Inertia matrix
I = np.array([[Ixx, -Ixy, -Ixz],
              [-Ixy, Iyy, -Iyz],
              [-Ixz, -Iyz, Izz]])
inv_I = np.linalg.inv(I)

def desired_attitude(t):
    # Desired attitude profile
    desired_roll = 0
    desired_pitch = -t * 2 * np.pi / (110 * 60)  # One full revolution over one orbit
    desired_yaw = 0
    return np.array([desired_roll, desired_pitch, desired_yaw])

def disturbance(t):
    Solar_pressure_torque = maximum_torque(q, A_s, c_ps, c_g)
    Gravity_gradient_torque = gravity_torque(r=r, theta=theta, I_y=Iyy, I_z=Izz)
    Aerodynamic_torque = aerodynamic_torque(C_d=C_d, rho=rho, A=A_aero, r=r, c_g=c_g, c_pa=c_pa)
    if t < 1000:
        Internal_disturbance_torque = np.array([0, 0, 0])
    # elif t >= 1000 and t < 2200:
    #     Internal_disturbance_torque = np.array([-1.8e-3, 0, 0]) # moving the antenna
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

def controller(omega, attitude_error, Kp, Kd):
    Kp = np.array([1.5, Kp, 1.5])   # Proportional gains for long term stability: [1, 1, 1]
    Kd = np.array([34.5, Kd, 34.5])      # Derivative gains for long term stability: [10, 15, 10]
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

def dynamics_factory(Kp, Kd):
    def dynamics(t, state):
        omega = state[0:3]
        attitude = state[3:6]
        attitude_error = attitude - desired_attitude(t)
        tau = controller(omega, attitude_error, Kp, Kd)
        d = disturbance(t)

        omega_cross = np.array([
            [0, -omega[2], omega[1]],
            [omega[2], 0, -omega[0]],
            [-omega[1], omega[0], 0]
        ])

        omega_dot = inv_I @ (tau + d - omega_cross @ (I @ omega))
        attitude_dot = omega  # Or your improved version
        return np.concatenate([omega_dot, attitude_dot])
    return dynamics

# Initial state: omega = [0,0,0], attitude = [0,0,0]
omega0 = np.array([om_x, om_y, om_z])
attitude0 = np.array([0, 0, 0])
state0 = np.concatenate([omega0, attitude0])

# Desired angular velocity
omega_desired = np.array([om_x, om_y, om_z])

# Time span and evaluation points
t_span = (0, 110 * 60) 
t_eval = np.linspace(t_span[0], t_span[1], 5000)

# Define kp and kd grids
kp_values = np.arange(0, 10, 0.5)         # 0 to 99
kd_values = np.arange(0, 50, 0.5)    # 0, 10, ..., 990

# Preallocate arrays to hold metrics
max_dev_roll = np.zeros((len(kp_values), len(kd_values)))
max_dev_pitch = np.zeros_like(max_dev_roll)
max_dev_yaw = np.zeros_like(max_dev_roll)
total_torque_x = np.zeros_like(max_dev_roll)
total_torque_y = np.zeros_like(max_dev_roll)
total_torque_z = np.zeros_like(max_dev_roll)
torque_activation = np.zeros_like(max_dev_roll)

for i, kp in enumerate(kp_values):
    for j, kd in enumerate(kd_values):
        Kp_vec = np.array([kp]*3)
        Kd_vec = np.array([kd]*3)
        dyn = dynamics_factory(kp, kd)

        sol = solve_ivp(dyn, t_span, state0, t_eval=t_eval, method='RK45')

        # Extract attitude angles
        roll = sol.y[3]
        pitch = sol.y[4]
        yaw = sol.y[5]

        # Desired attitudes for each timestep
        desired_att = np.array([desired_attitude(t) for t in sol.t]).T  # shape (3, len(t))

        # Index to ignore initial transient
        idx_start = 500

        # Max deviations in mrad
        max_dev_roll[i, j] = np.max(np.abs(roll[idx_start:] - desired_att[0, idx_start:])) * 1000
        max_dev_pitch[i, j] = np.max(np.abs(pitch[idx_start:] - desired_att[1, idx_start:])) * 1000
        max_dev_yaw[i, j] = np.max(np.abs(yaw[idx_start:] - desired_att[2, idx_start:])) * 1000

        # Calculate total torque used
        dt = sol.t[1] - sol.t[0]
        torque_sum_x = 0
        torque_sum_y = 0
        torque_sum_z = 0
        torque_activation_count = 0

        for k in range(len(sol.t)):
            omega_k = sol.y[0:3, k]
            attitude_k = sol.y[3:6, k]
            attitude_err_k = attitude_k - desired_att[:, k]
            tau_k = controller(omega_k, attitude_err_k, kp, kd)
            torque_sum_x += np.abs(tau_k[0]) * dt
            torque_sum_y += np.abs(tau_k[1]) * dt
            torque_sum_z += np.abs(tau_k[2]) * dt
            for torque in tau_k:
                if torque != 0:
                    torque_activation_count += 1

        total_torque_x[i, j] = torque_sum_x
        total_torque_y[i, j] = torque_sum_y
        total_torque_z[i, j] = torque_sum_z
        torque_activation[i, j] = torque_activation_count

        
       
    print(f'completed {(i+1)/len(kp_values)*100:.2f}% of Kp values')

# --- Prepare DataFrame for CSV export ---

Kp_grid, Kd_grid = np.meshgrid(kp_values, kd_values, indexing='ij')

data_dict = {
    "Kp": Kp_grid.flatten(),
    "Kd": Kd_grid.flatten(),
    "MaxDevRoll_mrad": max_dev_roll.flatten(),
    "MaxDevPitch_mrad": max_dev_pitch.flatten(),
    "MaxDevYaw_mrad": max_dev_yaw.flatten(),
    "TotalTorqueX_Nm_s": total_torque_x.flatten(),
    "TotalTorqueY_Nm_s": total_torque_y.flatten(),
    "TotalTorqueZ_Nm_s": total_torque_z.flatten(),
    "TorqueActivation": torque_activation.flatten()
}

df = pd.DataFrame(data_dict)

# Save CSV
df.to_csv("gain_sweep_results3.csv", index=False)
print("Saved gain sweep results to gain_sweep_results3.csv")