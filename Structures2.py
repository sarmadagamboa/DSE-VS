import math
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

class Material:
    def __init__(self, E, rho, s_yld, s_ult, alpha=0.8, n=0.6, nu=0.334):
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

class Polygon:
    def __init__(self, height, points, px, py, stringers, thickness):
        '''Initializes the Polygon class
        INPUTS:
            height: FLOAT -> height of the polygon
            points: LIST -> list of points defining the polygon
            px: LIST -> list of x-coordinates of the polygon points
            py: LIST -> list of y-coordinates of the polygon points
            stringers: LIST -> list of number of stringers per element
            thickness: FLOAT or NUMPY ARRAY -> thickness (range) of the polygon skin

        OUTPUTS:
            None
        '''
        self.height = height
        self.points = points
        self.px = px
        self.py = py
        self.stringers = stringers
        self.t = thickness

        self.element_lengths()
        
    def element_lengths(self):
        '''Calculates the elements of the polygon
        INPUTS:
            None
        OUTPUTS:
            None
        '''
        self.elements = []
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
        self.lengths = lengths
        self.material = material

        if self.type == 'hat':
            self.hat_area()
            self.hat_stress()
    
    def hat_stress(self):
        s_side_r = self.material.alpha * ((0.423/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[0])**2)**(1-self.material.n)
        s_vert_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[1])**2)**(1-self.material.n)
        s_top_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                         * (self.t/self.lengths[2])**2)**(1-self.material.n)

        if s_side_r > 1:
            s_side = self.material.s_yld
        else:
            s_side = s_side_r * self.material.s_yld
        
        if s_vert_r > 1:
            s_vert = self.material.s_yld
        else:
            s_vert = s_vert_r * self.material.s_yld
        
        if s_top_r > 1:
            s_top = self.material.s_yld
        else:
            s_top = s_top_r * self.material.s_yld
        
        a_side = self.t * self.lengths[0]
        a_vert = self.t * self.lengths[1]
        a_top = self.t * self.lengths[2]

        self.crip_stress = (2 * s_side * a_side + 2 * s_vert * a_vert + s_top * a_top) / (2 * a_side + 2 * a_vert + a_top)

    def hat_area(self):
        self.area = (2 * self.lengths[0] + 2 * self.lengths[1] + self.lengths[2]) * self.t + 4 * self.t**2


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

        
    def launcher_loading(self):
        fnat_ax = 20
        fnat_lat = 6
        maxg_ax = 6
        maxg_lat = 2

        p_ax = maxg_ax  * self.mass * 9.81
        p_lat =  maxg_lat * self.mass * 9.81

        p_eq = (p_ax + 2*(p_lat*self.geometry.height/2)/1.414) *1.25 #will be checked against both buckling and yield/ultimate strength
        ## 1.414 value above must be switched to half of the longest diagonal of the polygon
        
        #print(p_ax,p_lat,(p_lat*l/2),p_eq/1.25,p_eq)
        self.load = p_eq

        return self.load

    def calculate_loads(self):
        if isinstance(self.geometry, Polygon):
            # Calculate loads for polygon
            return self.calculate_poly_loads()
        elif isinstance(self.geometry, Cylinder):
            # Calculate loads for cylinder
            pass
        else:
            raise ValueError("Invalid geometry type")
    
    def calculate_poly_loads(self):
        self.w_e = (self.geometry.t/2) * math.sqrt((4.0 * math.pi**2)/(12.0 * (1 - self.material.nu**2))) * math.sqrt(self.material.E/self.stringer.crip_stress)
            
        print(f"w_e: {self.w_e}")

        for i in range(len(self.geometry.elements)):
            for n in range(len(self.geometry.t)):
                if self.geometry.elements[i] < 2 * self.w_e[n] * (self.geometry.stringers[i]):
                    print(f"Element {i} has too many stringers for thickness {self.geometry.t[n]}")
            
        self.element_empty_length = []
        for i in range(len(self.geometry.elements)):
            self.element_empty_length.append((self.geometry.elements[i] - 2 * self.w_e * (self.geometry.stringers[i]))/(self.geometry.stringers[i]+1))
            for n in range(len(self.geometry.t)):
                if self.element_empty_length[i][n] < 0:
                    self.element_empty_length[i][n] = 0

        self.element_empty_stress = []
        for i in range(len(self.geometry.elements)):
            self.element_empty_stress.append(4.0 * ((math.pi**2*self.material.E) / (12.0 * (1 - self.material.nu**2))) * (self.geometry.t / self.element_empty_length[i])**2)

            
        self.element_crippling_stress = []
        self.element_area = []
            
        for i in range(len(self.geometry.elements)):
            empty_area = self.element_empty_length[i] * (self.geometry.stringers[i]+1) * self.geometry.t
            empty_stress = self.element_empty_stress[i]
            stringered_area = self.geometry.elements[i] * self.geometry.t - empty_area
            stringer_area = self.stringer.area * self.geometry.stringers[i]
            stringer_stress = self.stringer.crip_stress
            self.element_crippling_stress.append((empty_area * empty_stress + (stringered_area + stringer_area) * stringer_stress)/(empty_area + stringered_area + stringer_area))
            self.element_area.append(empty_area + stringered_area + stringer_area)
                
        self.total_load_bearing = 0
        self.total_area = 0
        for i in range(len(self.geometry.elements)):
            self.total_load_bearing += self.element_crippling_stress[i] * self.element_area[i]
            self.total_area += self.element_area[i]
        self.total_stress = self.total_load_bearing / self.total_area

        #print(f"Total load bearing: {self.total_load_bearing}")
        #print(f"Total area: {self.total_area}")
        #print(f"Total stress: {self.total_stress}")
        print(self.total_load_bearing)
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




if __name__ == "__main__":
    points = [0, 1, 2, 3]
    x = [0, 2, 2, 0]
    y = [0, 0, 2, 2]

    stringers = [0, 0, 0, 0] #number of stringers per element
    #print(stringers)
    t = np.linspace(0.0001, 0.006, 50)
    #print(t)
    aluminium = Material(E=70e9, rho=2800, s_yld=448e6, s_ult=524e6) #based on aluminium 7075
    hat_stringer = Stringer(type='hat', thickness=0.0005, lengths=[0.01, 0.01, 0.01], material=aluminium)
    box = Polygon(height=1.2, points=points, px=x, py=y, stringers=stringers, thickness=t)
    load_calc = Load_calculation(geometry=box, stringer=hat_stringer, material=aluminium, sc_mass=691)

    #load_calc.launcher_loading()
    launchload = load_calc.launcher_loading()
    loadbearing = load_calc.calculate_loads()
    weight = load_calc.calculate_weight()

    print(f"Launch load: {launchload}")
    print(f"Load bearing: {loadbearing}")
    print(f"Weight: {weight}")


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