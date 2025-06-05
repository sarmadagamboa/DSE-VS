import numpy as np
from SSmodel import create_state_space_system
import control
import matplotlib.pyplot as plt

sys = create_state_space_system()

t = np.arange(0, 110*60, 1)  # Time vector for simulation
u = np.zeros((sys.B.shape[1], len(t)))  # Input matrix for control and disturbances

for i in range(len(t)):
    if i in range(5000, 5010):
        u[6, i] = 100 # strong "random" internal disturbance

u[3, len(t)-1] = 10 # constant aerodynamic torque

for i in range(len(t)):
    t_out, y = control.forced_response(sys, T=t, U=u[:, i:i+10])
    if y[0] > 0.01:
        u[0, i] = -1  # Apply control input in the x direction if angular velocity exceeds threshold
    elif y[0] < -0.01:
        u[0, i] = 1 # Apply control input in the x direction if angular velocity exceeds threshold


# control input in the x direction
u[0, :len(t)] = -1  # Example control input for the 

# control input in the y direction
u[1, :len(t)] = 0 

# control input in the z direction
u[2, :len(t)] = 0

# aerodynamic disturbance input
u[3, :len(t)] = 1  # Example disturbance input

# gravity gradient disturbance input
u[4, :len(t)] = 0  # Example disturbance input

# solar radiation pressure disturbance input
u[5, :len(t)] = 0  # Example disturbance input

# internal disturbance input
u[6, :len(t)] = 0  # Example disturbance input

t_out, y = control.forced_response(sys, T=t, U=u)

plt.figure()
plt.plot(t_out, y[0, :], label='Omega_x')
plt.plot(t_out, y[1, :], label='Omega_y')
plt.plot(t_out, y[2, :], label='Omega_z')
plt.title('Angular Velocity Response')
plt.xlabel('Time (s)')
plt.ylabel('Angular Velocity (rad/s)')
plt.legend()
plt.grid()
plt.show()