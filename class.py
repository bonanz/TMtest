# -*- coding: utf-8 -*-
"""
Spyder Editor
PartIII project
Class definition for a optical filter modeled using transfer matrix 
Layer thickness and wavelength need to be in the same unit. Usually nm is used
"""

import numpy as np


class layers:
    """
    define a object representing a multilayer filter
    """
    def __init__(self, N, d):
        """
        Initialise the object by supplying a list of refractive index and a list of thickness
        d in the unit the same as  the wavelength aka nms for convenience. \n
        N: a list containing the refractive indices \n
        d: a list containing the thickness of layers \n
        The films are indexed from the BACKSIDE
        """
        # store the stucture, note the refractive index can be complex
        self.N = np.array(N, dtype = complex)
        self.d = np.array(d, dtype = float)
    
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
        calcuate the overall trasfer matrices for each layer. 
        The result is stored as a stack of N matrix given by (N,M,M)
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
        # Store current wavelenght of potential reference
        self.Lambda = Lambda
    
    def calc_M_total(self):
        """
        Calculate the overall transfer matrix by taking dot products of matrix of
        each layer
        """
        self.M_total_s = self.M_s[:,:,0]
        self.M_total_p = self.M_p[:,:,0]
        for i in range(1,self.N.size):
            self.M_total_s = np.dot(self.M_total_s, self.M_s[:,:,i])
            self.M_total_p = np.dot(self.M_total_p, self.M_p[:,:,i])
            
    def calc_r(self, n_inc, n_ex, Lambda ,inc_ang = 0):
        """
        Calculate reflectivity of incident light with given angle
        n_inc: refractive index of the incident medium
        n_ex: refractive index of the exiting medium
        Lambda: wavelength in nm
        inc_ang: angle of incidence, the default is 0(normal incident)
        
        The function return a tuple of R,T. R and T are arrays containing two elements
        for two difference polarisation, e.g. R = [R_s, R_p].
        
        The function will exicute calculation of the transfer matrices of each layer
        and then take the dot producs of them.
        """
        Alpha = n_inc*np.sin(inc_ang/180*np.pi)
        self.calc_M(Alpha, Lambda)
        self.calc_M_total()
        M = np.dstack((self.M_total_s, self.M_total_p))
        # calculate pesodu indices for two sides
        eta_inc_s = np.sqrt(n_inc**2 - Alpha**2)
        eta_inc_p = n_inc**2 / eta_inc_s
        eta_ex_s = np.sqrt(n_ex**2 - Alpha**2)
        eta_ex_p = n_ex**2 / eta_ex_s
        # group indices together
        eta_inc = np.array([eta_inc_s, eta_inc_p])
        eta_ex = np.array([eta_ex_s, eta_ex_p])
        # Calculate r and t
        r = (eta_inc*M[0,0] - eta_ex*M[1,1] + eta_inc*eta_ex*M[0,1] - M[1,0]) \
          / (eta_inc*M[0,0] + eta_ex*M[1,1] + eta_inc*eta_ex*M[0,1] + M[1,0])
        t = 2*eta_inc / (eta_inc*M[0,0] + eta_ex*M[1,1] + eta_inc*eta_ex*M[0,1] + M[1,0])
        # Obtain power reflectivity and transmitivity
        R = np.absolute(r)**2
        T = np.absolute(t)**2*np.real(n_ex)/np.real(n_inc)
        return R,T

# Test if the class is working        
if __name__ == "__main__":
    R = []
    I = []
    # Test with calculation of the spectrum response
    for i in range(200,1000):
        l = layers([2,3,2],[50,100,50])
        R += [l.calc_r(1,1,i,60)[0]]
        I += [i]
    R = np.asarray(R)
    import matplotlib.pyplot as plt
<<<<<<< HEAD
    plt.plot(I,R[:,0], label = "s-polarisation" )
    plt.plot(I,R[:,1], label = "p-polarisation")
    plt.legend()
=======
    plt.plot(I,R)
    #new comment 2
>>>>>>> b06205bf1b732c2b23a9b2835d55bb68811dea85
