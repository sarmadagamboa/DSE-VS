def compute_MMOI(l_z = 1.2, w_y= 1.7, h_x = 3, sc_mass = 934): 
    I_z = 1/12 * sc_mass * (w_y**2 + h_x**2)
    I_x = 1/12 * sc_mass * (w_y**2 + l_z**2)
    I_y = 1/12 * sc_mass * (l_z**2 + h_x**2)

    I_xy = -1/24 * sc_mass * (w_y*h_x)
    I_xz = -1/24 * sc_mass * (l_z*h_x)
    I_yz = -1/24 * sc_mass * (w_y*l_z)

    return [I_x, I_y, I_z, I_xy, I_xz, I_yz]

# Call the function and print the result
mmoi_list = compute_MMOI()
print(mmoi_list)
