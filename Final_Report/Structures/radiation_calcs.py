def surface_area(length, width, height):
    return 2 * 2 * (length * width + length * height + width * height)

def shielding_mass(area, thickness_m, rho=900):
    volume = area * thickness_m  # m³
    return volume * rho  # kg

#rad_plastic = rad_si * 0.6 
Orbit_duration = 1.9 #in yrs 
#Mars_rad_Mgyday = 0.3 #mGy/day
Mars_rad_radyear_plastic = 10.95 #rad/year (plastic)
Mars_rad_radyear_si = Mars_rad_radyear_plastic/0.6 #krad/year (plastic)
print("plastic per year", Mars_rad_radyear_plastic, "Si (krad) per year", Mars_rad_radyear_si)
#commerical components: 5-10 rad plastic (grays) total 
#Mars_rad_radyear = 150 #rad/year (Si)
Mars_rad_radtotal = Orbit_duration*Mars_rad_radyear_plastic

Al_reduction = 35/100
Felt_rad_radtotal = (1-Al_reduction) * Mars_rad_radtotal
print("Total plastic per year, with Al reduction:", Felt_rad_radtotal)

limit_components = 10
reduction_required = (Felt_rad_radtotal - limit_components)  / (Felt_rad_radtotal) 
print(reduction_required)


#1 cm thickness HDPE required (rho = 900)


# Component dimensions (L, W, H) in meters
components = {
    "LRP": (0.25, 0.15, 0.25),
    "OBE": (0.3, 0.2, 0.07),
    #"PCDU": (0.007**(1/3),0.007**(1/3),0.007**(1/3))
}

# Parameters
thickness_cm = 1  # HDPE thickness in cm
thickness_m = thickness_cm / 100  # convert to meters
rho_hdpe = 900  # kg/m³

# Calculate and print results
print(f"HDPE Shielding Mass for {thickness_cm} cm thickness:\n")
rad_mass_tot = 0 
for name, dims in components.items():
    area = surface_area(*dims)
    mass = shielding_mass(area, thickness_m, rho_hdpe)
    print(f"{name}: Surface Area = {area:.4f} m², HDPE Mass = {mass:.3f} kg")
    rad_mass_tot += mass

print(rad_mass_tot)