# -*- coding: utf-8 -*-
"""
Spyder Editor
PartIII project
Class definition for a optical filter
This is a temporary script file.
"""

import numpy as np


class layers:
    """
    define a object representing a multilayer filter
    """
    def __init__(self, N, d):
        """
        Initialise the object by supplying a list of refractive index and a list of thickness
        d in the unit the same as  the wavelength aka nms for convenience.
        """
        # store the stucture
        self.N = np.array(N)
        self.d = np.array(d)
    
    def calc_ps(self, Alpha, Lambda):
        """
        calculate the phase shift term and store as a list of N elements
        """
        ps = 2*np.pi/Lambda*(np.sqrt(self.N**2-Alpha**2))*self.d
        self.ps = np.array(ps)
        
    def calc_pi(self, Alpha):
        """
        calculate the pesdo indices by the Alpha supplied
        pi[0]: list of s polarisation pesodu indices
        """
        pi = [0,0]
        # for s polarisation
        pi[0] = np.sqrt(self.N**2 - Alpha**2)
        # p polariastion
        pi[1] = self.N**2/pi[0]
        self.pi = np.array(pi)
        
    def calc_M(self, Alpha, Lambda):
        """
        calcuate the overall trasfer matrices for each layer
        the result is stored as a stack of N matrix given by (N,M,M)
        """
        self.calc_pi(Alpha)
        self.calc_ps(Alpha, Lambda)
        # start with s polarisation
        ms11 = np.cos(self.ps)
        ms12 = (np.complex(0,1)/self.pi[0])*np.sin(self.ps)
        ms21 = ms12*self.pi[0]**2
        ms22 = ms11
        self.M_s = np.array([[ms11, ms12], [ms21, ms22]], dtype = complex)
        # now do p polarisations
        mp12 = np.complex(0,1)/self.pi[1]*np.sin(self.ps)
        mp21 = mp12*self.pi[1]**2
        self.M_p = np.array([[ms11, mp12], [mp21, ms22]], dtype = complex)
    
    def calc_M_total(self):
        self.M_total_s = self.M_s[:,:,0]
        self.M_total_p = self.M_p[:,:,0]
        for i in range(1,self.N.size):
            self.M_total_s = np.dot(self.M_total_s, self.M_s[:,:,i])
            self.M_total_p = np.dot(self.M_total_p, self.M_p[:,:,i])
    def calc_r(self, n_inc, n_ex, Lambda ,inc_ang = 0):
        """
        Calculate reflectivity of incident light with given angle
        """
        Alpha = n_inc*np.sin(inc_ang/180*np.pi)
        self.calc_M(Alpha, Lambda)
        self.calc_M_total()
        M = np.dstack((self.M_total_s, self.M_total_p))
        # calculate pesodu indices
        eta_inc_s = np.sqrt(n_inc**2 - Alpha**2)
        eta_inc_p = n_inc**2 / eta_inc_s
        eta_ex_s = np.sqrt(n_ex**2 - Alpha**2)
        eta_ex_p = n_ex**2 / eta_ex_s
        eta_inc = np.array([eta_inc_s, eta_inc_p])
        eta_ex = np.array([eta_ex_s, eta_ex_p])
        r = (eta_inc*M[0,0] - eta_ex*M[1,1] + eta_inc*eta_ex*M[0,1] - M[1,0]) \
          / (eta_inc*M[0,0] + eta_ex*M[1,1] + eta_inc*eta_ex*M[0,1] + M[1,0])
        t = 2*eta_inc / (eta_inc*M[0,0] + eta_ex*M[1,1] + eta_inc*eta_ex*M[0,1] + M[1,0])
        R = np.absolute(r)**2
        T = np.absolute(t)**2*np.real(n_ex)/np.real(n_inc)
        return R,T
if __name__ == "__main__":
    R = []
    I = []
    for i in range(200,1000):
        l = layers([1,2,3,5],[70,90,100,120])
        R += [l.calc_r(1,1,i,0)[1][0]]
        I += [i]
    import matplotlib.pyplot as plt
    plt.plot(I,R)