import math
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

class Material: #add other properties such as cost or manufacturability 
    def __init__(self, name, E, rho, s_yld, s_ult, cost_kg, manuf, alpha, n, nu, alpha_therm, min_realistic_t):
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
        self.name = name 
        self.E = E
        self.rho = rho
        self.s_yld = s_yld
        self.s_ult = s_ult
        self.alpha = alpha
        self.cost_kg = cost_kg
        self.manuf = manuf 
        self.n = n
        self.nu = nu
        self.alpha_therm = alpha_therm
        self.min_realistic_t = min_realistic_t
        

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
        self.width = self.px[1]-self.px[0]
        self.length = self.py[2]-self.py[1]
        #self.halfmaxdiag = 0.5 * 

        self.element_lengths()
        self.max_diagonal()
        
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

    def max_diagonal(self):
        # maximum diagonal (non-edge) of the polygon
        max_diag = 0 
        for i in range(len(self.points)):
            for j in range(i+1, len(self.points)):
                # Skip adjacent points (they form edges)
                if abs(i - j) == 1 or (i == 0 and j == len(self.points) - 1):
                    continue
                dx = self.px[i] - self.px[j]
                dy = self.py[i] - self.py[j]
                diag = math.sqrt(dx**2 + dy**2)
                if diag > max_diag:
                    max_diag = diag
        self.max_diag = max_diag


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
    def __init__(self, type, thickness, lengths, material, manuf_stringer):
        self.type = type
        self.t = thickness
        self.lengths = lengths #lengths of each of the sides 
        self.material = material
        self.manuf_stringer = manuf_stringer 

        if self.type == 'hat':
            self.hat_area()
            self.hat_stress()
        elif self.type == 'Z':
            self.z_area()
            self.z_stress()
        elif self.type == 'I':
            self.i_area()
            self.i_stress()

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
    
    def hat_radius_of_gyration(self):
        t = self.t
        l0, l1, l2 = self.lengths  # side flange, vertical web, top plate
        
        # Segment areas
        a_side = t * l0
        a_vert = t * l1
        a_top = t * l2
        a_corner = t * t  # corner reinforcement patches
        
        # y-positions of centroids from base (assuming flat base)
        y_side = t / 2                     # base flange center
        y_vert = l1 / 2 + t               # web center
        y_top = l1 + t + t / 2            # top plate center
        y_corner_bottom = t / 2              # corners (if placed just above base)
        y_corner_top = t + l1 + t/2

        # Area-weighted centroid (neutral axis)
        A_total = 2 * a_side + 2 * a_vert + a_top + 4 * a_corner
        y_bar = (2 * a_side * y_side + 2 * a_vert * y_vert + a_top * y_top + 2 * a_corner * y_corner_bottom + 2 * a_corner * y_corner_top) / A_total

        # got here - Moments of inertia about the centroidal axis
        I_total = 0
        # Sides (base flanges)
        I_total += 2 * ((1/12) * l0 * t**3 + a_side * (y_bar - y_side)**2)
        # Verticals
        I_total += 2 * ((1/12) * t * l1**3 + a_vert * (y_bar - y_vert)**2)
        # Top plate
        I_total += (1/12) * l2 * t**3 + a_top * (y_bar - y_top)**2
        # Corners
        I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_top)**2)
        I_total += 2 * ((1/12) * t * t**3 + a_corner * (y_bar - y_corner_bottom)**2)

        # Radius of gyration
        r = (I_total / A_total) ** 0.5
        self.r_gyr = r
    
    def z_stress(self):
        s_side_r = self.material.alpha * ((0.423/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[0])**2)**(1-self.material.n)
        s_vert_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[1])**2)**(1-self.material.n)
        s_top_r = self.material.alpha * ((0.423/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[0])**2)**(1-self.material.n)

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

        self.crip_stress = (s_side * a_side + s_vert * a_vert + s_top * a_top) / (a_side + a_vert + a_top)

    def z_area(self):
        self.area = (self.lengths[0] + self.lengths[1] + self.lengths[2]) * self.t + 2 * self.t**2

    def i_stress(self):
        s_side_r = self.material.alpha * ((0.423/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[0])**2)**(1-self.material.n)
        s_vert_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[1])**2)**(1-self.material.n)

        if s_side_r > 1: #yielding happens first 
            s_side = self.material.s_yld
        else: #crippling happes first, s_side is then either the yielding or critical depending on which is most constraining 
            s_side = s_side_r * self.material.s_yld
        
        if s_vert_r > 1:
            s_vert = self.material.s_yld
        else:
            s_vert = s_vert_r * self.material.s_yld
        

        a_side = self.t * self.lengths[0] #not truly necessary 
        a_vert = self.t * self.lengths[1]

        self.crip_stress = (4 * s_side * a_side + s_vert * a_vert) / (4* a_side + a_vert)

    def i_area(self):
        self.area = (4 * self.lengths[0] + self.lengths[1]) * self.t + 2 * self.t**2

class Load_calculation:
    def __init__(self, geometry, stringer, material, sc_mass,t):
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
        self.t = t

        
    def launcher_loading(self, fnat_ax=20, fnat_lat=6, maxg_ax = 6, maxg_lat = 2):
        self.fnat_ax = fnat_ax
        self.fnat_lat = fnat_lat
        self.maxg_ax = maxg_ax
        self.maxg_lat = maxg_lat

        p_ax = maxg_ax  * self.mass * 9.81
        p_lat =  maxg_lat * self.mass * 9.81


        # CHANGE 
        p_eq = (p_ax + 2*(p_lat*self.geometry.height/2)/self.geometry.max_diag) *1.25 #will be checked against both buckling and yield/ultimate strength
        ## TODO!!! 1.414 value above must be switched to half of the longest diagonal of the polygon
        
        
        #print(p_ax,p_lat,(p_lat*l/2),p_eq/1.25,p_eq)
        self.peq_load = p_eq
        print(p_eq)
        return self.peq_load

    def calculate_loads(self,kd_factor):
        if isinstance(self.geometry, Polygon):
            # Calculate loads for polygon
            return self.calculate_poly_loads(kd_factor)
        elif isinstance(self.geometry, Cylinder):
            # Calculate loads for cylinder
            pass
        else:
            raise ValueError("Invalid geometry type")
    
    def calculate_poly_loads(self, kd_factor): #stringer crippling stress is the average 
        #w_e is actually an array, calculting a different w_e value for each different thickness value within the thickness array 
        #which is why there is a loop over n as well for each of the thicknee-dependent values 
        self.w_e = (self.geometry.t/2) * math.sqrt((4* math.pi**2)/(12.0 * (1 - self.material.nu**2))) * math.sqrt(self.material.E/self.stringer.crip_stress)
        print("Crippling stress:", self.stringer.crip_stress)   
        #print(f"w_e: {self.w_e}")

        """ 
        for i in range(len(self.geometry.elements)):
            for n in range(len(self.geometry.t)):
                
                #if self.geometry.elements[i] < 2 * self.w_e[n] * (self.geometry.stringers[i]): #stringers array has n. of stringers per edge, so for each edge we check 
                    #print(f"Element {i} has too many stringers for thickness {self.geometry.t[n]}")
                #print(f"w_e: {self.w_e}")
                if self.geometry.elements[i] < 2 * self.w_e[n] * (self.geometry.stringers[i]):
                    #print("ERROR:", self.geometry.t[n], self.geometry.elements[i], 2 * self.w_e[n] * (self.geometry.stringers[i]))
                    print(f"i={i}, material = {self.material.name}, thickness = {self.geometry.t[n]} element length={self.geometry.elements[i]:.4f}, stringers={self.geometry.stringers[i]}, w_e={self.w_e[n]:.4f}, required length={2*self.w_e[n]*self.geometry.stringers[i]:.4f}")
                    raise ValueError(f"Invalid configuration: Element {i} has too many stringers ({self.geometry.stringers[i]}) and we ({self.w_e[n]}) for thickness {self.geometry.t[n]:.6f} m")
        """



        self.element_empty_length = [] #calculates for each element, the empty parts without stringers, i.e., not just the total for each element, but the smaller parts within each element, i.e., unstiffened region between each stringer
        for i in range(len(self.geometry.elements)):
            self.element_empty_length.append((self.geometry.elements[i] - 2 * self.w_e * (self.geometry.stringers[i]))/(self.geometry.stringers[i]+1))
            for n in range(len(self.geometry.t)):
                if self.element_empty_length[i][n] < 0: #stringers 
                    self.element_empty_length[i][n] = 0

        self.element_empty_stress = [] #calculates the stress of the empty skin (for each element) 
        for i in range(len(self.geometry.elements)):
            self.element_empty_stress.append(kd_factor * 4.0 * ((math.pi**2*self.material.E) / (12.0 * (1 - self.material.nu**2))) * (self.geometry.t / self.element_empty_length[i])**2)

            
        self.element_crippling_stress = []
        self.element_area = []
            
        for i in range(len(self.geometry.elements)): #for each element. each of these are lists depending on the t  
            empty_area = self.element_empty_length[i] * (self.geometry.stringers[i]+1) * self.geometry.t
            empty_stress = self.element_empty_stress[i]
            stringered_area = self.geometry.elements[i] * self.geometry.t - empty_area #area on panel considering the w_e
            stringer_area = self.stringer.area * self.geometry.stringers[i] #actual area of the hat stringers 
            stringer_stress = self.stringer.crip_stress * kd_factor #each stringer crippling stress based on the previous weighted average of each of the sides 
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


    def calculate_weight(self,safety_mass_margin):
        '''Calculates the weight of the structure
        INPUTS:
            None
        OUTPUTS:
            None
        '''
        
        self.weight = self.total_area * self.material.rho * self.geometry.height
        self.weight *= safety_mass_margin
        return self.weight

    def check_frequency(self, min_t_crippling, min_t_ix):
        A_req_tot = (self.fnat_ax/0.25)**2*self.mass*self.geometry.height/self.material.E
        I_req_tot = (self.fnat_lat/0.56)**2*self.mass*self.geometry.height**3/self.material.E

        #print(self.geometry.length, self.geometry.width)

        I_skin = min_t_crippling*(self.geometry.length**2*self.geometry.width + self.geometry.width**2*self.geometry.length)*1/3
        I_req_stringers = I_req_tot - I_skin
        A_factor_list = [] 
        d = [] 
        #print(self.geometry.elements)
        for i in range(len(self.geometry.elements)): #for each element. each of these are lists depending on the t  
            if i == 0 or i == 2: 
                Ad2_A_factor = self.geometry.stringers[i] * (self.geometry.length/2)**2
                A_factor_list.append(Ad2_A_factor) 
            if i == 1 or i == 3: 
                stringers_half = self.geometry.stringers[i]/2
                length_between = self.geometry.length/(stringers_half+1)

                for j in range(1, int(stringers_half) + 1):
                    #print(j)
                    d = length_between * j 
                    Ad2_A_factor = d**2
                    Ad2_A_factor = d**2 * 2 #for both sides
                    A_factor_list.append(Ad2_A_factor) 
                
            #sum all factors (elements) in the A_factor list 
        Ad2_A_factor_sum = sum(A_factor_list)
        required_stringer_area =  I_req_stringers/Ad2_A_factor_sum
        
        
        if required_stringer_area > self.stringer.area: #only one value here, same for all thicknesses
            raise ValueError("Stringer area insufficient for inertial/bending requirement") 
        else: 
            print("Sufficient stringer area")
        #needs to be checked for minimum case out of all thicknesses w
        # otherwise self.total_area varies for all thicknesses within each stringer combination/loop, would be checked for all thicknesses for each call of this
        
        if A_req_tot > self.total_area[min_t_ix]: 
            raise ValueError("Total area insufficient for area requirement")
        else: 
            print("Sufficient total area")

        """
        t_Areq = A_req/(2*self.geometry.elements[0] + 2*self.geometry.elements[1])
        t_Ireq = 3*I_req/(self.geometry.elements[0]**2*self.geometry.elements[1] + self.geometry.elements[1]**2*self.geometry.elements[0])

        t_sreq_ult = self.peq_load / (2 * (self.geometry.elements[0] + self.geometry.elements[1]) * self.material.s_ult) * 1.1
        t_sreq_yld = self.peq_load / (2 * (self.geometry.elements[0] + self.geometry.elements[1]) * self.material.s_yld) * 1.1

        self.rigidity_t = max(t_Areq, t_Ireq, t_sreq_ult, t_sreq_yld)
        print(t_Areq, t_Ireq, t_sreq_ult, t_sreq_yld, self.rigidity_t, self.geometry.elements[0], self.geometry.elements[1])
        """
         

    def calculate_t_min(self, loadbearing, launchload): 
        satisfying_t_indices_cr = np.where(loadbearing >= launchload)[0] #0 purely bc of syntax
        #print(satisfying_t_indices_cr)
        # minimum crippling thickness
        if len(satisfying_t_indices_cr) == 0:
            min_t_crippling = None
            #print("No thickness satisfies the load bearing requirement.")
            
        else:
            min_t_crippling = t[satisfying_t_indices_cr[0]] #gets the minimum by accessing the first index of the possible ts
            min_t_ix = satisfying_t_indices_cr[0]
        
        return min_t_crippling, min_t_ix


    def check_thermal(self, instrument_path, delta_T_max, max_shift):
        delta_L = self.material.alpha_therm * instrument_path * delta_T_max
        valid_delta_L = (delta_L <= max_shift)

        #s_therm = self.material.E * self.material.alpha_therm * delta_T_max / ((1-self.material.nu)) 
        
        if valid_delta_L == False: 
            raise ValueError("Thermal requirements not met")
        else: 
            print("Thermal requirements met")
        #conservative, without division by 2, to account for stringers' effect  
        #s_therm = self.material.E * self.material.alpha_therm * delta_T_max / ((1-self.material.nu))
    
    def calculate_johnson(self): 
        s_johnson_crippling = self.stringer.stringer_stress(1 - (self.stringer.stringer_stress * (self.geometry.height*self.stringer.r_gyr)**2)/(4 * math.pi()**2 * self.material.E) )
        euler_buckling_yielding = math.pi()**2 * self.material.E * (self.stringer.r_gyr/self.geometry.height)**2
        critical_slenderness_ratio = math.sqrt(2 * math.pi()**2 * self.material.E /self.stringer.stringer_stress) 
        #then based on that calculate limit 
    

if __name__ == "__main__":
    ### INPUTS ###
    x = [0, 1.7, 1.7, 0] #x-coordinates of the polygon points
    y = [0, 0, 1.2, 1.2] #y-coordinates of the polygon points
    #stringers = [0, 0, 0, 0] #number of stringers per element
    sc_mass = 1063 #mass of the total wet spacecraft in kg
    sc_height = 3 #height of the spacecraft in m
    

    #introducing realism, change  
    kd_factor = 0.85 #knockdown factor 
    safety_mass_margin = 1.15 #for joints, discontinuities etc.

    #thermal stability inputs 
    instrument_path = sc_height
    delta_T_max = 10
    max_shift = 1.5e-3

    best_weight = float('inf')
    best_config = None
    w_mass = 0.7          # mass is critical for launch
    w_cost = 0.3          # cost still matters
    w_manuf = 0.1        # manufacturability is important, but less so than the other two

    ### CREATE OBJECTS ###
    #aluminium = Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    
    materials = [
    Material(name="Aluminium 7075", E=70e9, rho=2800, s_yld=448e6, s_ult=524e6, cost_kg=5.0, manuf=0.9, alpha=0.8, n=0.6, nu=0.334, alpha_therm = 23.6e-6, min_realistic_t = 0.0015),  # Aluminium 7075
    ]
    #Material(name="Titanium Ti-6Al-4V", E=116e9, rho=4420, s_yld=880e6, s_ult=950e6, cost_kg=30.0, manuf=0.7, alpha=0.75, n=0.15, nu=0.34, alpha_therm = 8.6e-6, min_realistic_t = 0.0020),  # Titanium Ti-6Al-4V
    #Material(name="Stainless Steel", E=230e9, rho=7850, s_yld=250e6, s_ult=460e6, cost_kg=1.5, manuf=0.5, alpha=0.7, n=0.25, nu=0.3, alpha_therm = 17.2e-6, min_realistic_t = 0.0030),   # Stainless Steel
    #Material(name="CFRP", E=150e9, rho=1600, s_yld=800e6, s_ult=1100e6, cost_kg=90.0, manuf=0.4, alpha=0.6, n=0.05, nu=0.2, alpha_therm = 0.5e-6, min_realistic_t = 0.0015),  # CFRP
    


    possible_configs = [] 
    stringer_count_fail_configs = [] 
    frequency_fail_configs = [] 
    thermal_fail_configs = [] 
    
    stringers_iterations_4 = [[2, 2, 2, 2]]                                                    


    for material in materials: 
        
        t = np.linspace(0.0015, 0.0023, 100)

        stringer_types = [
        #Stringer(type='hat', thickness=0.005, lengths=[0.03,0.03,0.03], material=material, manuf_stringer = 0.7),
        Stringer(type='Z', thickness=0.0028, lengths=[0.05,0.05,0.05], material=material, manuf_stringer = 0.9),
        #Stringer(type='I', thickness=0.006, lengths=[0.03,0.03], material=material, manuf_stringer = 1),
        ] 


        for stringer in stringer_types: 
            stringer_type_iter = stringer
    
            for stringers_iter in stringers_iterations_4:
                    box = Polygon(height=sc_height, px=x, py=y, stringers=stringers_iter, thickness=t)
                    load_calc = Load_calculation(geometry=box, stringer=stringer_type_iter, material=material, sc_mass=sc_mass, t=t)
                    
                    # rest of your processing (launcher loading, load-bearing, weight, etc.)

                    launchload = load_calc.launcher_loading() #string 
                    
                    """
                    ###### crippling check, determination of loadbearing criterion, checks for too many stringers per panel in calculation 
                    try: 
                        loadbearing = load_calc.calculate_loads(kd_factor = kd_factor)
                    except ValueError as e:
                        stringer_count_fail_configs.append({
                        "material": material,
                        "material name": material.name,
                        "stringers": stringers_iter,
                        "stringer type": stringer.type,
                        })
                        print(f"Too many stringers per panel. Skipping invalid config {stringers_iter}: {e}")
                        print()
                        continue
                    weight = load_calc.calculate_weight(safety_mass_margin = safety_mass_margin)
                    min_t_crippling, min_t_ix= load_calc.calculate_t_min(loadbearing, launchload)
                    """
                    loadbearing = load_calc.calculate_loads(kd_factor = kd_factor)
                    weight = load_calc.calculate_weight(safety_mass_margin = safety_mass_margin)
                    min_t_crippling, min_t_ix= load_calc.calculate_t_min(loadbearing, launchload)
                    ###### frequency check 
                    try:
                        
                        load_calc.check_frequency(min_t_crippling, min_t_ix)
                    except ValueError as e:
                        frequency_fail_configs.append({
                        "material": material,
                        "material name": material.name,
                        "stringers": stringers_iter,
                        "stringer type": stringer.type,
                        })
                        print(f"Insufficient dimensions for this configuration. Skipping invalid config {stringers_iter}: {e}")
                        continue 
                    
                    ###### thermal check 
                    try: 
                        load_calc.check_thermal(instrument_path, delta_T_max, max_shift)
                    except ValueError as e:
                        thermal_fail_configs.append({
                        "material": material,
                        "material name": material.name,
                        "stringers": stringers_iter,
                        "stringer type": stringer.type,
                        })
                        print(f"Thermal expansion/stresses too high. Skipping invalid config {stringers_iter}: {e}")
                        continue 
                    
                    
                    ##### finalise calculations
                    #weight_t_min = weight[min_t_ix]
                    weight_t_min = weight[min_t_ix] + 0.75 + 8.370000000000001 + 10.56 + 3.2082468 #ADCS plate + radiation shielding + MMOD protection + CAI vibration
                    cost = weight_t_min * material.cost_kg

                    # Normalize individual metrics
                    inv_mass_score = 1 / weight_t_min
                    inv_cost_score = 1 / cost
                    manuf_score = material.manuf * stringer_type_iter.manuf_stringer

                    # Weighted combined score
                    #score = (w_mass * inv_mass_score) + (w_cost * inv_cost_score) + (w_manuf * manuf_score)
                    score = 1/weight_t_min

                    possible_configs.append({
                        "material": material,
                        "material name": material.name,
                        "stringers": stringers_iter,
                        "stringer type": stringer.type,
                        "t_min": min_t_crippling,
                        "mass": weight_t_min,
                        "cost": cost,
                        "score": score,
                        "launch load to bear": launchload,
                    })

                    """
                    if weight_t_min < best_weight:
                        best_weight = weight_t_min
                        best_config = (stringers_iter, min_t_crippling, best_weight)
                    """
                    """ 
                    print(f"Configuration: Stringers per edge: {stringers_iter}") 
                    print(f"Launch load: {launchload}")
                    print(f"Load bearing: {loadbearing}")
                    print(f"Minimum thickness (for load-bearing): {min_t_crippling:.6f} m")
                    print(f"Minimum corresponding mass: {weight_t_min:.6f} kg")
                    print()
                    possible_configs.append(best_config)
                    """
                    """

            if best_config:
                best_stringers_iter, best_min_t_crippling, best_best_weight = best_config
                print(f"Best configuration:\n Stringers per edge: {best_stringers_iter}\n Minimum thickness: {best_min_t_crippling:.6f} m\n Minimum mass: {best_best_weight:.2f} kg")
            else:
                print("No feasible configuration found.")
        """
    print(possible_configs)

    # Sort and print best configs
    possible_configs.sort(key=lambda r: r["score"], reverse=True) 

    print("Top material configurations for Mars mission:")
    for r in possible_configs[:5]: #get first top configs 
        m = r["material"]
        print(f"Material Name={m.name}, E={m.E/1e9:.0f} GPa, Cost=${m.cost_kg}/kg, Manufacturability={m.manuf}")
        print(f"Stringers: {r['stringers']}, Stringer type: {r['stringer type']}, t_min: {r['t_min']:.4f} m, Mass: {r['mass']:.2f} kg, Cost: ${r['cost']:.2f}, Score: {r['score']:.2e}\n")
        total_mass =r['mass']

    print("Stringer fail configs:", stringer_count_fail_configs) 
    print("Frequency fail configs:", frequency_fail_configs) 
    print("Thermal fail configs:", thermal_fail_configs) 

    #Final mass 

    print("Final mass:", total_mass)
    