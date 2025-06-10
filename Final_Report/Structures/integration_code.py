from Structure2_optimal_current2 import * 
from compute_MMOI import * 

structural_mass, dim_height, dim_length, dim_width, sc_mass = structural_mass_wpanel(sc_mass = 921.51, dim_height = 3, dim_length = 1.2, dim_width = 1.7)
structural_mass_itercode = structural_mass_wpanel(sc_mass = 921.51, dim_height = 3, dim_length = 1.2, dim_width = 1.7)[0]
MMOI_list = compute_MMOI(l_z = dim_length, w_y= dim_width, h_x = dim_height, sc_mass = sc_mass)

print(structural_mass, MMOI_list, structural_mass_itercode)