
m_prop = 365
ox_fuel_ratio = 1.65
final_press = 2000000 #Pa
storage_temp = 300#K
storage_press = 20000000 #Pa
gas_const = 2077 # J/kg*K
def compute_mass_fuel(m_prop, ox_fuel_ratio):
    fuel_mass = m_prop/(1+ox_fuel_ratio)
    oxidiser_mass =ox_fuel_ratio*fuel_mass
    return fuel_mass, oxidiser_mass
def compute_pressurant(final_press,storage_press, storage_temp, gas_const, v_prop_biprop):
    M_press = final_press * v_prop_biprop/(gas_const*storage_temp)
    v_press = M_press*gas_const*storage_temp/storage_press
    return M_press, v_press

fuel_mass, oxidiser_mass = compute_mass_fuel(m_prop,ox_fuel_ratio)
print(fuel_mass, oxidiser_mass)
v_fuel_required = fuel_mass / 880
v_oxidizer_required = oxidiser_mass / 1370
v_prop_biprop = v_fuel_required + v_oxidizer_required
compute_pressurant(final_press,storage_press, storage_temp, gas_const, v_prop_biprop)
pressurant_mass,pressurant_volume = compute_pressurant(final_press,storage_press, storage_temp, gas_const, v_prop_biprop)
print(pressurant_mass)