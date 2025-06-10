#lateral resonance critical 
#A_req_tot = (self.fnat_ax/0.25)**2*self.mass*self.geometry.height/self.material.E4
br_height = 0.01
br_width = 0.493
br_length = 0.493
I_req = (100/0.56)**2*br_height**3/(4*10e9)
   
#t_Areq = A_req/(2*br_width + 2*br_length)
I_actual = (br_length * br_width**3) / 12
print(I_req, I_actual)
V = br_height * br_width * br_length 
mass = V * 1320 
print(mass)
#t_Ireq = 3*I_req/(br_width**2*br_length + br_length**2*br_width)

#t_sreq_ult = self.peq_load / (2 * (self.geometry.elements[0] + self.geometry.elements[1]) * self.material.s_ult) * 1.1
#t_sreq_yld = self.peq_load / (2 * (self.geometry.elements[0] + self.geometry.elements[1]) * self.material.s_yld) * 1.1

#self.rigidity_t = max(t_Areq, t_Ireq, t_sreq_ult, t_sreq_yld)
#print(t_Areq, t_Ireq, t_sreq_ult, t_sreq_yld, self.rigidity_t, self.geometry.elements[0], self.geometry.elements[1])
        