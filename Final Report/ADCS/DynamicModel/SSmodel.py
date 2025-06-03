import numpy as np
from ModelConstants import *
import control

def Omega():
    # Angular velocity vector
    # This is the value that the system will try to maintain
    omega = np.array([[om_x],
                     [om_y],
                     [om_z]])
    return omega

def Omega_dot():
    # Angular acceleration vector
    omega_dot = np.array([[0],
                         [0],
                         [0]])
    return omega_dot

def Omega_cross(omega):
    # skew-symmetric matrix for the angular velocity vector
    omega_cross = np.array([[0, -omega[2][0], omega[1][0]],
                           [omega[2][0], 0, -omega[0][0]],
                           [-omega[1][0], omega[0][0], 0]])
    return omega_cross

def I_mat():
    # Inertia matrix
    I = np.array([[Ixx, -Ixy, -Ixz],
              [-Ixy, Iyy, -Iyz],
              [-Ixz, -Iyz, Izz]])
    return I

def Tau():
    # Control torque vector. First column is the control torque in the x direction, 
    # second column is the control torque in the y direction, 
    # and third column is the control torque in the z direction
    # The numbers in this matrix is the minimum torque bits for the control torque
    tau = np.array([[min_x_bit, 0, 0],
                    [0, min_y_bit, 0],
                    [0, 0, min_z_bit]])
    return tau

def Dist():
    # disturbance torque vector. First column is the aerodynamic torque, 
    # second column is the gravity gradient torque, 
    # and third column is the solar radiation pressure torque, 
    # and the fourth column is internal disturbances
    # The numbers in this matrix are the constant parts of the disturbance torque, 
    # their time varying parts are included in the input vector
    dist = np.array([[ae_x, grav_x, sol_x, int_x],
                  [ae_y, grav_y, sol_y, int_y],
                  [ae_z, grav_z, sol_z, int_z]])
    return dist

def create_state_space_system():
    # Create the state space system
    omega = Omega()
    omega_dot = Omega_dot()
    omega_cross = Omega_cross(omega)
    I_matrix = I_mat()
    tau = Tau()
    dist = Dist()

    # Dynamic model of the spacecraft
    C1 = I_matrix
    C2 = np.dot(omega_cross, I_matrix)
    C3 = np.hstack((tau, dist))

    # State space representation of the dynamic model
    A = -np.dot(np.linalg.inv(C1), C2)
    B = -np.dot(np.linalg.inv(C1), C3)
    C = np.eye(A.shape[0])
    D = np.zeros((A.shape[0], B.shape[1]))

    # Create the state space system
    sys = control.ss(A, B, C, D)
    
    return sys
