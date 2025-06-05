import math
from scipy.optimize import fsolve
from itertools import product 
import numpy as np
import matplotlib.pyplot as plt

class Material: #add other properties such as cost or manufacturability 
    def __init__(self, E, rho, s_yld, s_ult, alpha, n, nu, cost, manuf):
        '''Initializes the Material class
        INPUTS:
            E: FLOAT -> Young's modulus in Pa
            rho: FLOAT -> density in kg/m^3
            s_yld: FLOAT -> yield stress in Pa
            s_ult: FLOAT -> ultimate stress in Pa
            alpha, n, nu are correct for aluminium alloys
        
        OUTPUTS:
            None
        '''
        self.E = E
        self.rho = rho
        self.s_yld = s_yld
        self.s_ult = s_ult
        self.alpha = alpha
        self.n = n
        self.nu = nu
        self.cost = cost 
        self.manuf = manuf

class Polygon: #polygon-shaped cross-section of structural member
    def __init__(self, height, px, py, stringers, thickness):
        '''Initializes the Polygon class
        INPUTS:
            height: FLOAT -> height of the polygon
            points: LIST -> list of points defining the polygon
            px: LIST -> list of x-coordinates of the polygon points
            py: LIST -> list of y-coordinates of the polygon points
            stringers: LIST -> list of number of stringers per element (per edge)
            thickness: FLOAT or NUMPY ARRAY -> thickness (range) of the polygon skin, to be able to test different thicknesses

        OUTPUTS:
            None
        '''
        self.height = height
        self.points = list(range(0, len(px)))
        self.px = px
        self.py = py
        self.stringers = stringers
        self.t = thickness

        self.element_lengths()
        
    def element_lengths(self):
        '''Calculates the length of each polygon side through Pythagoras (assuming straight polygon edges) 
        INPUTS:
            None
        OUTPUTS:
            None
        '''
        self.elements = [] #element lengths 
        for i in range(len(self.points)):
            next_i = (i + 1) % len(self.points)
            length = math.sqrt((self.px[i] - self.px[next_i])**2 + (self.py[i] - self.py[next_i])**2)
            self.elements.append(length)

class Cylinder:
    '''Initializes the Cylinder class
        INPUTS:
            height: FLOAT -> height of the cylinder
            radius: FLOAT -> radius of the cylinder
            stringers: FLOAT -> number of stringers on the cylinder
            thickness: FLOAT -> thickness of the cylinder skin

        OUTPUTS:
            None
        '''
    def __init__(self, height, radius, stringers, thickness):
        self.height = height
        self.radius = radius
        self.stringers = stringers
        self.t = thickness



class Stringer:
    def __init__(self, type, thickness, lengths, material):
        self.type = type
        self.t = thickness
        self.lengths = lengths #lengths of each of the sides 
        self.material = material

        if self.type == 'hat':
            self.hat_area()
            self.hat_stress()
        elif self.type == 'Z':
            self.z_area()
            self.z_stress()
        #add other types 
    
    def hat_stress(self): #fraction formulae here 
        s_side_r = self.material.alpha * ((0.423/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[0])**2)**(1-self.material.n)
        s_vert_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[1])**2)**(1-self.material.n)
        s_top_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                         * (self.t/self.lengths[2])**2)**(1-self.material.n)

        if s_side_r > 1: #yielding happens first 
            s_side = self.material.s_yld
        else: #crippling happes first, s_side is then either the yielding or critical depending on which is most constraining 
            s_side = s_side_r * self.material.s_yld
        
        if s_vert_r > 1:
            s_vert = self.material.s_yld
        else:
            s_vert = s_vert_r * self.material.s_yld
        
        if s_top_r > 1:
            s_top = self.material.s_yld
        else:
            s_top = s_top_r * self.material.s_yld
        
        a_side = self.t * self.lengths[0] #not truly necessary 
        a_vert = self.t * self.lengths[1]
        a_top = self.t * self.lengths[2]
        
        #average crippling stress based on area-weighted average
        self.crip_stress = (2 * s_side * a_side + 2 * s_vert * a_vert + s_top * a_top) / (2 * a_side + 2 * a_vert + a_top)

    def hat_area(self): #sides and corners 
        self.area = (2 * self.lengths[0] + 2 * self.lengths[1] + self.lengths[2]) * self.t + 4 * self.t**2

    def z_stress(self):
        t, l0, l1 = self.t, *self.lengths
        E, nu, alpha, n, s_yld = self.material.E, self.material.nu, self.material.alpha, self.material.n, self.material.s_yld
        def crippling(l):
            return alpha * ((4/E * (math.pi**2 * E) / (12 * (1 - nu**2)) * (t/l)**2)**(1-n)) * s_yld
        s_flange = min(s_yld, crippling(l0))
        s_web = min(s_yld, crippling(l1))
        a_flange = t * l0 * 2
        a_web = t * l1
        self.crip_stress = (s_flange * a_flange + s_web * a_web) / (a_flange + a_web)

    def z_area(self):
        self.area = (2 * self.lengths[0] + self.lengths[1]) * self.t + 2 * self.t**2


class Load_calculation:
    def __init__(self, geometry, stringer, material, sc_mass):
        '''Initializes the Load_calculation class
        INPUTS:
            geometry: OBJECT -> geometry of the structure
            stringer: OBJECT -> stringer properties
            material: OBJECT -> material properties

        OUTPUTS:
            None
        '''
        self.geometry = geometry
        self.stringer = stringer
        self.material = material
        self.mass = sc_mass

        
    def launcher_loading(self, fnat_ax=20, fnat_lat=6, maxg_ax = 6, maxg_lat = 2):
        self.fnat_ax = fnat_ax
        self.fnat_lat = fnat_lat
        self.maxg_ax = maxg_ax
        self.maxg_lat = maxg_lat

        p_ax = maxg_ax  * self.mass * 9.81
        p_lat =  maxg_lat * self.mass * 9.81

        # CHANGE 
        p_eq = (p_ax + 2*(p_lat*self.geometry.height/2)/1.7614) *1.25 #will be checked against both buckling and yield/ultimate strength
        ## TODO!!! 1.414 value above must be switched to half of the longest diagonal of the polygon
        
        #print(p_ax,p_lat,(p_lat*l/2),p_eq/1.25,p_eq)
        self.peq_load = p_eq

        return self.peq_load

    def calculate_loads(self):
        if isinstance(self.geometry, Polygon):
            # Calculate loads for polygon
            return self.calculate_poly_loads()
        elif isinstance(self.geometry, Cylinder):
            # Calculate loads for cylinder
            pass
        else:
            raise ValueError("Invalid geometry type")
    
    def calculate_poly_loads(self): #stringer crippling stress is the average 
        #w_e is actually an array, calculting a different w_e value for each different thickness value within the thickness array 
        #which is why there is a loop over n as well for each of the thicknee-dependent values 
        self.w_e = (self.geometry.t/2) * math.sqrt((4.0 * math.pi**2)/(12.0 * (1 - self.material.nu**2))) * math.sqrt(self.material.E/self.stringer.crip_stress)
            
        print(f"w_e: {self.w_e}")

        for i in range(len(self.geometry.elements)):
            for n in range(len(self.geometry.t)):
                if self.geometry.elements[i] < 2 * self.w_e[n] * (self.geometry.stringers[i]): #stringers array has n. of stringers per edge, so for each edge we check 
                    print(f"Element {i} has too many stringers for thickness {self.geometry.t[n]}")
            
        self.element_empty_length = [] #calculates for each element, the empty parts without stringers, i.e., not just the total for each element, but the smaller parts within each element, i.e., unstiffened region between each stringer
        for i in range(len(self.geometry.elements)):
            self.element_empty_length.append((self.geometry.elements[i] - 2 * self.w_e * (self.geometry.stringers[i]))/(self.geometry.stringers[i]+1))
            for n in range(len(self.geometry.t)):
                if self.element_empty_length[i][n] < 0: #stringers 
                    self.element_empty_length[i][n] = 0

        self.element_empty_stress = [] #calculates the stress of the empty skin (for each element) 
        for i in range(len(self.geometry.elements)):
            self.element_empty_stress.append(4.0 * ((math.pi**2*self.material.E) / (12.0 * (1 - self.material.nu**2))) * (self.geometry.t / self.element_empty_length[i])**2)

            
        self.element_crippling_stress = []
        self.element_area = []
            
        for i in range(len(self.geometry.elements)): #for each element. each of these are lists depending on the t  
            empty_area = self.element_empty_length[i] * (self.geometry.stringers[i]+1) * self.geometry.t
            empty_stress = self.element_empty_stress[i]
            stringered_area = self.geometry.elements[i] * self.geometry.t - empty_area #area on panel considering the w_e
            stringer_area = self.stringer.area * self.geometry.stringers[i] #actual area of the hat stringers 
            stringer_stress = self.stringer.crip_stress #each stringer crippling stress based on the previous weighted average of each of the sides 
            #overall crippling stress weighted sum calculated per element 
            self.element_crippling_stress.append((empty_area * empty_stress + (stringered_area + stringer_area) * stringer_stress)/(empty_area + stringered_area + stringer_area))
            self.element_area.append(empty_area + stringered_area + stringer_area)
                
        self.total_load_bearing = 0
        self.total_area = 0
        for i in range(len(self.geometry.elements)):
            self.total_load_bearing += self.element_crippling_stress[i] * self.element_area[i] #actual load P, which is why only multiplied by the A
            self.total_area += self.element_area[i]
        self.total_stress = self.total_load_bearing / self.total_area #weighted average of all of the elements for entire structure 

        #print(f"Total load bearing: {self.total_load_bearing}")
        #print(f"Total area: {self.total_area}")
        #print(f"Total stress: {self.total_stress}")
        return self.total_load_bearing

    def calculate_weight(self):
        '''Calculates the weight of the structure
        INPUTS:
            None
        OUTPUTS:
            None
        '''
        self.weight = self.total_area * self.material.rho

        print(f"Weight: {self.weight}")

        return self.weight

    def calculate_area_inertia_thickness(self):
        A_req = (self.fnat_ax/0.25)**2*self.mass*self.geometry.height/self.material.E
        I_req= (self.fnat_lat/0.56)**2*self.mass*self.geometry.height**3/self.material.E
        t_Areq = A_req/(2*self.geometry.elements[0] + 2*self.geometry.elements[1])
        t_Ireq = 3*I_req/(self.geometry.elements[0]**2*self.geometry.elements[1] + self.geometry.elements[1]**2*self.geometry.elements[0])

        t_sreq_ult = self.peq_load / (2 * (self.geometry.elements[0] + self.geometry.elements[1]) * self.material.s_ult) * 1.1
        t_sreq_yld = self.peq_load / (2 * (self.geometry.elements[0] + self.geometry.elements[1]) * self.material.s_yld) * 1.1

        self.rigidity_t = max(t_Areq, t_Ireq, t_sreq_ult, t_sreq_yld)
        print(t_Areq, t_Ireq, t_sreq_ult, t_sreq_yld, self.rigidity_t, self.geometry.elements[0], self.geometry.elements[1])

        return self.rigidity_t

    

    def evaluate_structure(self):
        weight = self.calculate_weight()
        margin = self.stringer.crip_stress / self.peq_load
        return {'weight': weight, 'margin': margin, 'valid': margin >= 1.0}


if __name__ == "__main__":
    ### INPUTS ###
    x = [0, 1.45, 1.45, 0] #x-coordinates of the polygon points
    y = [0, 0, 1, 1] #y-coordinates of the polygon points
    stringers = [0, 0, 0, 0] #number of stringers per element
    t = np.linspace(0.0001, 0.006, 50) #thickness (range) of the skin
    sc_mass = 1063 #mass of the total wet spacecraft in kg
    sc_height = 4.5 #height of the spacecraft in m
    str_type = ['hat']

    ###### CREATE OBJECTS #####
    # MATERIAL ITERATION for now just iterating through materials 
    aluminium = Material("Al7075", 70e9, 2800, 448e6, 524e6, cost=1.0, manufacturability=0.9)
    titanium = Material("Ti6Al4V", 113.8e9, 4430, 880e6, 950e6, cost=1.8, manufacturability=0.6)
    cfrp = Material("CFRP", 130e9, 1600, 600e6, 800e6, alpha=0.6, n=0.55, nu=0.3, cost=2.2, manufacturability=0.7)
    materials = [aluminium, titanium, cfrp]

    # STRINGER ITERATION , will then iterate through stringer properties 
    hat_stringer = Stringer(type=str_type, thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)

    # SHAPE ITERATION 
    box = Polygon(height=sc_height, px=x, py=y, stringers=stringers, thickness=t)

    load_calc = Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass)

    #load_calc.launcher_loading()
    launchload = load_calc.launcher_loading()
    loadbearing = load_calc.calculate_loads()

    weight = load_calc.calculate_weight()
    rigidity_t = load_calc.calculate_area_inertia_thickness()

    print(f"Launch load: {launchload}")
    print(f"Load bearing: {loadbearing}")
    print(f"Weight: {weight}")
    print(f"Rigidity thickness: {rigidity_t}")

    # Plot setup
    fig, ax1 = plt.subplots()

    # Plot load bearing on the left y-axis
    ax1.plot(t, loadbearing, 'b-', label='Load Bearing')
    ax1.set_xlabel('Thickness (t)')
    ax1.set_ylabel('Load Bearing (N)', color='b')
    ax1.tick_params(axis='y', labelcolor='b')

    ax1.axhline(y=launchload, color='g', linestyle='--', label=f'Limit: {launchload} N')

    # Create a second y-axis for mass
    ax2 = ax1.twinx()
    ax2.plot(t, weight, 'r-', label='Mass')
    ax2.set_ylabel('Mass (kg)', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # Add legends and grid
    fig.tight_layout()
    plt.title('Load Bearing and Mass vs. Thickness')
    plt.grid(True)
    plt.show()

    # Automatically find minimum thickness 
    # Find minimum thickness satisfying load bearing ≥ launch load
    satisfying_t_indices = np.where(loadbearing >= launchload)[0] #0 purely bc of syntax

    if len(satisfying_t_indices) == 0:
        min_t_bearing = None
        print("No thickness satisfies the load bearing requirement.")
    else:
        min_t_bearing = t[satisfying_t_indices[0]] #gets the minimum by accessing the first index of the possible ts
        weight_min_t_bearing = weight[satisfying_t_indices[0]]

        min_t_bearing_overall = max(min_t_bearing, rigidity_t)
        print(f"Minimum thickness to satisfy load-bearing: {min_t_bearing_overall:.6f} m")
        
        print(f"Minimum corresponding mass FOR LOAD BEARING ONLY (not rigidity): {weight_min_t_bearing:.6f} kg")
    
    
    best_config = None
    best_score = float('inf')

    #for "iteration names" in "actual code array names"
    for (mat_skin, t_skin, typ) in product(materials, t, str_type): 
        stringer = Stringer(type=typ, thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=mat_skin)
        box = Polygon(height=sc_height, px=x, py=y, stringers=stringers, thickness=t_skin)
        load_calc = Load_calculation(geometry=box, stringer=stringer, material=mat_skin, sc_mass=sc_mass)
        result = load_calc.evaluate_structure()

        score = result['weight'] * mat_skin.manufacturability + mat_str.cost
        if result['valid'] and score < best_score:
            best_score = score
            best_config = (mat_skin.name, mat_str.name, t_s, t_str, typ, stringer.area, result['margin'], result['weight'])

    print(f"Best configuration:\n Skin Material: {best_config[0]}\n Stringer Material: {best_config[1]}\n Skin Thickness: {best_config[2]:.6f} m\n Stringer Thickness: {best_config[3]:.6f} m\n Type: {best_config[4]}\n Stringer Area: {best_config[5]:.6e} m^2\n Margin: {best_config[6]:.2f}\n Total Mass: {best_config[7]:.2f} kg")
    

    """ 
    stringer_thicknesses = [0.0003, 0.0005]
    stringer_counts = [[2, 2, 2, 2]]
    types = ['hat', 'Z']

    stringer_lengths = {
        'hat': [[0.01, 0.01, 0.01]],
        'Z': [[0.015, 0.01]]
    }

    best_config = None
    best_score = float('inf')

    for (mat_skin, mat_str, t_s, t_str, typ) in product(materials, materials, t_skin, stringer_thicknesses, types):
        for sl in stringer_lengths[typ]:
            for sc in stringer_counts:
                stringer = Stringer(type=typ, thickness=t_str, lengths=sl, material=mat_str)
                box = Polygon(height=sc_height, px=x, py=y, stringers=sc, thickness=[t_s]*4)
                calc = Load_calculation(geometry=box, stringer=stringer, material=mat_skin, sc_mass=sc_mass)
                result = calc.evaluate_structure()

                score = result['weight'] * mat_skin.manufacturability + mat_str.cost
                if result['valid'] and score < best_score:
                    best_score = score
                    best_config = (mat_skin.name, mat_str.name, t_s, t_str, typ, stringer.area, result['margin'], result['weight'])

    print(f"Best configuration:\n Skin Material: {best_config[0]}\n Stringer Material: {best_config[1]}\n Skin Thickness: {best_config[2]:.6f} m\n Stringer Thickness: {best_config[3]:.6f} m\n Type: {best_config[4]}\n Stringer Area: {best_config[5]:.6e} m^2\n Margin: {best_config[6]:.2f}\n Total Mass: {best_config[7]:.2f} kg")
    """