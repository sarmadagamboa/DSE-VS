import math
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

class Material: #add other properties such as cost or manufacturability 
    def __init__(self, E, rho, s_yld, s_ult, alpha=0.8, n=0.6, nu=0.334, alpha_therm = 23.6e-6):
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
        self.alpha_therm = alpha_therm

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
        p_eq = (p_ax + 2*(p_lat*self.geometry.height/2)/1.7614) *1.25 #will be checked against both buckling and yield/ultimate strength
        ## TODO!!! 1.414 value above must be switched to half of the longest diagonal of the polygon
        
        #print(p_ax,p_lat,(p_lat*l/2),p_eq/1.25,p_eq)
        self.peq_load = p_eq

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
        self.w_e = (self.geometry.t/2) * math.sqrt((4.0 * math.pi**2)/(12.0 * (1 - self.material.nu**2))) * math.sqrt(self.material.E/self.stringer.crip_stress)
            
        #print(f"w_e: {self.w_e}")

        for i in range(len(self.geometry.elements)):
            for n in range(len(self.geometry.t)):
                
                #if self.geometry.elements[i] < 2 * self.w_e[n] * (self.geometry.stringers[i]): #stringers array has n. of stringers per edge, so for each edge we check 
                    #print(f"Element {i} has too many stringers for thickness {self.geometry.t[n]}")
               
                if self.geometry.elements[i] < 2 * self.w_e[n] * (self.geometry.stringers[i]):
                    raise ValueError(f"Invalid configuration: Element {i} has too many stringers ({self.geometry.stringers[i]}) for thickness {self.geometry.t[n]:.6f} m")
 



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
        
        self.weight = self.total_area * self.material.rho
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


if __name__ == "__main__":
    ### INPUTS ###
    x = [0, 1.45, 1.45, 0] #x-coordinates of the polygon points
    y = [0, 0, 1, 1] #y-coordinates of the polygon points
    #stringers = [0, 0, 0, 0] #number of stringers per element
    min_realistic_t = 0.0025
    t = np.linspace(min_realistic_t, 0.006, 50) #thickness (range) of the skin
    sc_mass = 1063 #mass of the total wet spacecraft in kg
    sc_height = 4.5 #height of the spacecraft in m
    limiting_cr = 1

    #introducing realism, change  
    kd_factor = 0.85 #knockdown factor 
    safety_mass_margin = 2.5 #for joints, discontinuities etc.

    #thermal stability inputs 
    instrument_path = sc_height
    delta_T_max = 10
    max_shift = 1.5e-3

    best_weight = float('inf')
    best_config = None

    ### CREATE OBJECTS ###
    aluminium = Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = Stringer(type='hat', thickness=0.0005, lengths=[0.01,0.01,0.01], material=aluminium)
    
    possible_configs = [] 
    for s0 in range(0, 9, 2):
        for s1 in range(0, 9, 2):
            s2 = s0  # enforce symmetry
            s3 = s1
            stringers_iter = [s0, s1, s2, s3]

            box = Polygon(height=sc_height, px=x, py=y, stringers=stringers_iter, thickness=t)
            load_calc = Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass, t=t)
            
            launchload = load_calc.launcher_loading() #string 
            
            ###### crippling check, determination of loadbearing criterion, checks for too many stringers per panel in calculation 
            try: 
                loadbearing = load_calc.calculate_loads(kd_factor = kd_factor)
            except ValueError as e:
                print(f"Too many stringers per panel. Skipping invalid config {stringers_iter}: {e}")
                print()
                continue
            weight = load_calc.calculate_weight(safety_mass_margin = safety_mass_margin)
            min_t_crippling, min_t_ix= load_calc.calculate_t_min(loadbearing, launchload)

            ###### frequency check 
            try:
                load_calc.check_frequency(min_t_crippling, min_t_ix)
            except ValueError as e:
                print(f"Insufficient dimensions for this configuration. Skipping invalid config {stringers_iter}: {e}")
                continue 
            
            ###### thermal check 
            try: 
                load_calc.check_thermal(instrument_path, delta_T_max, max_shift)
            except ValueError as e:
                print(f"Thermal expansion/stresses too high. Skipping invalid config {stringers_iter}: {e}")
                continue 
            
            
            ##### finalise calculations
            weight_t_min = weight[min_t_ix]

            if weight_t_min < best_weight:
                best_weight = weight_t_min
                best_config = (stringers_iter, min_t_crippling, best_weight)
            
            
            print(f"Configuration: Stringers per edge: {stringers_iter}") 
            print(f"Launch load: {launchload}")
            print(f"Load bearing: {loadbearing}")
            print(f"Minimum thickness (for load-bearing): {min_t_crippling:.6f} m")
            print(f"Minimum corresponding mass: {weight_t_min:.6f} kg")
            print()
            possible_configs.append(best_config)


    if best_config:
        best_stringers_iter, best_min_t_crippling, best_best_weight = best_config
        print(f"Best configuration:\n Stringers per edge: {best_stringers_iter}\n Minimum thickness: {best_min_t_crippling:.6f} m\n Minimum mass: {best_best_weight:.2f} kg")
    else:
        print("No feasible configuration found.")

    
   

    """ 
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
    
    

    if limiting_cr == 0: #rigidity/frequency consideration is limiting 
                box = Polygon(height=sc_height, px=x, py=y, stringers=stringers_iter, thickness=np.array([min_t_overall]))
                load_calc = Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=sc_mass, t=t)

                # Recalculate load-bearing (needed to trigger internal area calculation)
                launchload = load_calc.launcher_loading() #string 
                try:
                    loadbearing = load_calc.calculate_loads(kd_factor = kd_factor)
                except ValueError as e:
                    print(f"Skipping invalid config {stringers_iter}: {e}")
                    print()
                    continue
                weight_t_min = load_calc.calculate_weight(safety_mass_margin =  safety_mass_margin)[0] 
    """