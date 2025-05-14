import math
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

class Material:
    def __init__(self, E, rho, f_yld, f_ult, alpha=0.8, n=0.6, nu=0.33):
        '''Initializes the Material class
        INPUTS:
            E: FLOAT -> Young's modulus in Pa
            rho: FLOAT -> density in kg/m^3
            s_yld: FLOAT -> yield stress in Pa
            s_ult: FLOAT -> ultimate stress in Pa
        
        OUTPUTS:
            None
        '''
        self.E = E
        self.rho = rho
        self.s_yld = f_yld
        self.s_ult = f_ult
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
            len = math.sqrt((self.px[i] - self.px[next_i])**2 + (self.py[i] - self.py[next_i])**2)
            self.elements.append(len)

class Cylinder:
    '''Initializes the Polygon class
        INPUTS:
            height: FLOAT -> height of the cylinder
            radius: FLOAT -> radius of the cylinder
            stringers: FLOAT -> number of stringers on the cylinder
            thickness: FLOAT -> thickness of the polygon skin

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
            self.hat_stress()
    
    def hat_stress(self):
        s_side_r = self.material.alpha * ((0.423/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[0])**2)**(1-self.material.n)
        s_vert_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                          * (self.t/self.lengths[1])**2)**(1-self.material.n)
        s_top_r = self.material.alpha * ((4/self.material.s_yld) * (math.pi**2 * self.material.E)/(12 * (1- self.material.nu**2))
                                         * (self.t/self.lengths[2])**2)**(1-self.material.n)

        if s_side_r > self.material.s_yld:
            s_side = self.material.s_yld
        else:
            s_side = s_side_r * self.material.s_yld
        
        if s_vert_r > self.material.s_yld:
            s_vert = self.material.s_yld
        else:
            s_vert = s_vert_r * self.material.s_yld
        
        if s_top_r > self.material.s_yld:
            s_top = self.material.s_yld
        else:
            s_top = s_top_r * self.material.s_yld
        
        a_side = self.t * self.lengths[0]
        a_vert = self.t * self.lengths[1]
        a_top = self.t * self.lengths[2]

        self.crip_stress = (2 * s_side * a_side + 2 * s_vert * a_vert + s_top * a_top) / (2 * a_side + 2 * a_vert + a_top)

        #### ADD AREA CALCULATION FOR THE STRINGER


class Load_calculation:
    def __init__(self, geometry, stringer, material):
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

        
    def launcher_loading(self): #####INCOMPLETE
        fnat_ax = 20
        fnat_lat = 6
        maxg_ax = 6
        maxg_lat = 2

        p_ax = maxg_ax  * m * 9.81
        p_lat =  maxg_lat * m * 9.81

        p_eq = (p_ax + 2*(p_lat*l/2)/r) *1.25 #will be checked against both buckling and yield/ultimate strength
        #print(p_ax,p_lat,(p_lat*l/2),p_eq/1.25,p_eq)
        self.load = p_eq

    def calculate_loads(self):
        if isinstance(self.geometry, Polygon):
            # Calculate loads for polygon
            self.calculate_poly_loads()
        elif isinstance(self.geometry, Cylinder):
            # Calculate loads for cylinder
            pass
        else:
            raise ValueError("Invalid geometry type")
    
    def calculate_poly_loads(self):
            w_e = (self.geometry.t/2) * math.sqrt((4.0 * math.pi**2)/(12.0 * (1 - self.material.nu**2))) * math.sqrt(self.material.E/self.stringer.crip_stress)

            for i in range(len(self.geometry.elements)):
                if self.geometry.elements[i] < 2* w_e * (self.geometry.stringers[i]-1):
                    print("Element {i} has too many stringers")
            
            self.element_empty_length = []
            for i in range(len(self.geometry.elements)):
                self.element_empty_length.append((self.geometry.elements[i] - 2 * w_e * (self.geometry.stringers[i]))/(self.geometry.stringers[i]+1))

            self.element_empty_stress = []
            for i in range(len(self.geometry.elements)):
                self.element_empty_stress.append(4.0*((math.pi**2)/(12.0 * 1 - self.material.nu**2)) * (self.geometry.t/self.element_empty_length[i])**2)
            






if __name__ == "__main__":
    aluminium = Material(70e9, 2800, 448e6, 524e6) #based on aluminium 7075

    points = [0, 1, 2, 3]
    x = [0, 3, 3, 0]
    y = [0, 0, 3, 3]

    stringers = []