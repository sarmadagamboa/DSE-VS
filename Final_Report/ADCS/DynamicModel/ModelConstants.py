import numpy as np


# Moments of inertia for the spacecraft
Ixx = 337
Iyy = 812.6
Izz = 925.4
Ixy = 0 #-289
Ixz = 0 #-199.3
Iyz = 0 #-64.2

# Initial angular velocity components in rad/s
om_x = 0
om_y = -2*np.pi/(110*60)  # one full revolution over one orbit
om_z = 0

# Minimum torque bits for the control torque in Nm
min_x_bit = 0.00005 # roughly the bits achieved by thrusters on GRACE
min_y_bit = 0.00005
min_z_bit = 0.00005

min_torque = -0.08
max_torque = 0.08


##### Disturbance torque constants #####
# Solar radiation pressure torque parameters
Antenna_diameter = 1.64
Antenna_depth = 0.3
SC_length = 3
SC_height = 1.2
SC_width = 1.7
Solar_array_area = 7.2
Solar_length = Solar_array_area/2/SC_length
worst_case_solar_incidence = 27/180*np.pi  # Worst case solar incidence angle in radians, assuming the solar panels are perpendicular to the sun
boom_length = np.tan(worst_case_solar_incidence)*SC_width
mounting_location = 1.2
q = 0.6 # preliminary reflectance factor, a more comprehensive analysis is needed as it will vary across the spacecraft
A_sp = Solar_array_area/2 # per solar array
A_side = SC_length*SC_height # Area of the side of the spacecraft
A_a = (Antenna_diameter/2)**2 * np.pi  # Area of the antenna
A_s = A_sp*2 + A_side
c_ps = (A_sp*Solar_length/2+A_sp*(Solar_length+2*boom_length+SC_height + Solar_length/2)+A_side*(Solar_length+boom_length+SC_height/2))/ A_s  # Location of center of solar pressure, Weighted average position of the center of solar pressure, assuming the antenna is perpendicular to the sun
c_g_s = Solar_length + boom_length + SC_height/2  # Location of center of gravity, Assuming the center of gravity is at the center of the spacecraft side, measured from the bottom of the solar panels
c_g = np.array([0, 0, 0]) # center of gravity is at the center of the spacecraft body, which is by deifinition at 0, 0, 0
c_ps = np.array([0, 0, c_g_s-c_ps])
# Aerodynamic torque parameters
C_d = 2.6 # Drag coefficient
rho = 1.4e-11  # Density of the Martian atmosphere in kg/m^3
r = 212.48e3  # Distance from the center of Mars to the spacecraft in meters
Solar_thickness = 10/1000
Sidepiece_height = 0.5
Sidepiece_width = 0.4
c_g_a = SC_width/2 # measured from the location where the arrays are mounted
A_aero = SC_height*SC_width + Antenna_depth*Antenna_diameter/2 + 2*Solar_length*Solar_thickness + Sidepiece_width*Sidepiece_height
c_pa = (SC_height*SC_width*SC_width/2 + Antenna_diameter*Antenna_depth/2*SC_width + Sidepiece_height*Sidepiece_width*(-Sidepiece_width/2)+2*Solar_length*Solar_thickness)/(A_aero)
c_pa = np.array([0, c_pa-c_g_a, 0.003])

print(f"C_pa = {c_pa} m and c_g_a = {c_g_a} m, measured from the bottom of the spacecraft")
# gravity gradient torque parameters
theta = 0 # Maximum torque occurs at pi/4, but this will never occur for the spacecraft