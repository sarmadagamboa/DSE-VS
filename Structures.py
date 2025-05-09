import numpy as np

# material properties, based on aluminium 2024-T6
E = 70e9  # Young's modulus in Pa
nu = 0.33  # Poisson's ratio
G = E / (2 * (1 + nu))  # shear modulus
rho = 7850  # density in kg/m^3	



# geometry
L = 4.6  # length in m
r = 1.5  # radius in m
b = 0.02  # stiffener spacing

# spacecraft properties
m = 1400  # mass in kg
fnat_axial = 20 # natural frequency in Hz
fnat_lateral = 6 # natural frequency in Hz


A_axial = ((fnat_axial/0.25)**2 * m * L)/E  # area needed axial
I_lateral = ((fnat_lateral/0.56)**2 * m * L**3)/E  # I needed lateral

print("Axial area needed: ", A_axial)
print("Lateral moment of inertia needed: ", I_lateral)

import numpy as np

# material properties, based on aluminium 2024-T6
E = 70e9  # Young's modulus in Pa
nu = 0.33  # Poisson's ratio
G = E / (2 * (1 + nu))  # shear modulus
rho = 7850  # density in kg/m^3	



# geometry
L = 4.6  # length in m
r = 1.5  # radius in m
b = 0.02  # stiffener spacing

# spacecraft properties
m = 1400  # mass in kg
fnat_axial = 20 # natural frequency in Hz
fnat_lateral = 6 # natural frequency in Hz


A_axial = ((fnat_axial/0.25)**2 * m * L)/E  # area needed axial
I_lateral = ((fnat_lateral/0.56)**2 * m * L**3)/E  # I needed lateral

print("Axial area needed: ", A_axial)
print("Lateral moment of inertia needed: ", I_lateral)

# material properties, based on aluminium 2024-T6
E = 70e9  # Young's modulus in Pa
nu = 0.33  # Poisson's ratio
G = E / (2 * (1 + nu))  # shear modulus
rho = 7850  # density in kg/m^3	



# geometry
L = 4.6  # length in m
r = 1.5  # radius in m
b = 0.02  # stiffener spacing

# spacecraft properties
m = 1400  # mass in kg
fnat_axial = 20 # natural frequency in Hz
fnat_lateral = 6 # natural frequency in Hz


A_axial = ((fnat_axial/0.25)**2 * m * L)/E  # area needed axial
I_lateral = ((fnat_lateral/0.56)**2 * m * L**3)/E  # I needed lateral

print("Axial area needed: ", A_axial)
print("Lateral moment of inertia needed: ", I_lateral)

# skin thickness based on axial-area requirement
t_axial = A_axial / (2 * np.pi * r)

# skin thickness based on bending-inertia requirement
t_bending = 

print(f"Thickness from axial requirement: {t_axial} m")
print(f"Thickness from bending requirement: {t_bending} m")

# choose the larger of the two
t_min = max(t_axial, t_bending)
print(f"Minimum required skin thickness: {t_min} m")
