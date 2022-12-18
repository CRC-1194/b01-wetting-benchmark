import numpy as np
import math
import ast

class Funcs:
    ''' ODE model based on
    https://doi.org/10.1016/j.jcp.2019.109221
    '''

    #function to solve ODE for contact angle evolution; returns the contact angle array
    @classmethod
    def get_ref_contact_angle(self, t):
        #**** VELOCITY MODEL FROM Mathis ADVECTION PAPER****#
        fx = lambda v0, t, tau, x1, x2: v0*(math.cos(np.pi*t / tau)) * ((-math.sin(np.pi*x1)) * (math.cos(np.pi*x2))) # ODE
        fy = lambda v0, t, tau, x1, x2: v0*(math.cos(np.pi*t / tau)) * ((math.cos(np.pi*x1)) * (math.sin(np.pi*x2))) # ODE

        #gradient matrix of the velocity [x0 y0: x1 y1]
        del_fx0 = lambda v0, t, tau, x1, x2: np.pi*v0*(math.cos(np.pi*t / tau)) * ((-math.cos(np.pi*x1))*(math.cos(np.pi*x2)))
        del_fx1 = lambda v0, t, tau, x1, x2: np.pi*v0*(math.cos(np.pi*t / tau)) * ((-math.sin(np.pi*x1))*(math.sin(np.pi*x2)))
        del_fy0 = lambda v0, t, tau, x1, x2: np.pi*v0*(math.cos(np.pi*t / tau)) * ((math.sin(np.pi*x1))*(math.sin(np.pi*x2)))
        del_fy1 = lambda v0, t, tau, x1, x2: np.pi*v0*(math.cos(np.pi*t / tau)) * ((math.cos(np.pi*x1))*(math.cos(np.pi*x2)))

        h = t[1]-t[0]

        v0 = 0.1
        t0 = 0.0 # start time
        tau = 0.2 # tau
        height0 = 0.1 #height at t=0
        radius0 = 0.2 #foot radius at t=0
        theta0 = math.acos(height0 / radius0) # CA at t=0

        #interface unnormalized normal at t=0
        lst_norm = [[-radius0*math.sin(theta0)],
                    [radius0*math.cos(theta0)]]
        interface_normal0 = np.array(lst_norm) 

        #intial position at the contact line at t=0
        lst_x0 = [[0.2], [0]]
        x0 = np.array(lst_x0)

        # Explicit Euler Method
        x = np.zeros((len(t),2))
        x[0 , 0] = x0[0]
        x[0 , 1] = x0[1]
    
        #interface_normal_not_normalized
        n =np.zeros((len(t),2))
        n[0 , 0] = interface_normal0[0]
        n[0 , 1] = interface_normal0[1]


        for i in range(0, len(t) - 1):
            x[i + 1, 0] = x[i ,0] + h*fx(v0, t[i], tau, x[i,0], x[i,1])
            x[i + 1, 1] = x[i ,1] + h*fy(v0, t[i], tau, x[i,0], x[i,1])


        for i in range(0, len(t) - 1):
            del_v = np.array([[del_fx0(v0, t[i], tau, x[i,0], x[i,1]), del_fy0(v0, t[i], tau, x[i,0], x[i,1])], 
                              [del_fx1(v0, t[i], tau, x[i,0], x[i,1]), del_fy1(v0, t[i], tau, x[i,0], x[i,1])]])
            del_v.transpose()
            n[i + 1] = n[i] - h * (del_v.dot(n[i]))


        #normalized normal
        n_hat =np.zeros((len(t),2))
        for i in range(0, len(n)):
            n_hat[i] = n[i] / np.linalg.norm(n[i])


        #boundary outward unit normal
        boundary_normal = np.array([[0], [-1]])

        contact_angle =np.zeros((len(t))) #ODE solution
        for i in range(0, len(contact_angle)):
            contact_angle[i] = math.acos(-np.dot(n_hat[i],boundary_normal)) *180 / np.pi

        return contact_angle


    @classmethod
    #Function to calculate the maximum error in contact angle
    def cal_max_error(self, ca_sim, ca_ref):
        diff = abs(ca_sim[0]-ca_ref[0])
        diff_temp = abs(ca_sim[0]-ca_ref[0])

        for i in range(len(ca_sim)):
            diff_temp = abs(ca_sim[i]-ca_ref[i])

            if (diff_temp > diff):
                diff = diff_temp

        return abs(diff)

    @classmethod
    #Function gives dictionary mapping variation number to parameter vector
    def fileMapping(self):
        var_map = {}
        var_lines = open("variation_file").readlines()
        nX =[]

        for line in var_lines:
        # Skip lines without mapping
            if '{' not in line:
                continue
            var_num = int(line.split()[1])
            dict_start = line.find('{')
            # Mappings in variation file can directly be interpreted by Python
            var_map[var_num] = ast.literal_eval(line[dict_start:-1])

        for key, value in var_map.items():
            meshSize = str(value).split(': ')[1][:-1]
            nX.append( int(1.0/float(meshSize)) )

        return nX
