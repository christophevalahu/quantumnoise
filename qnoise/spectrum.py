# QUANTUMNOISE
# Toolbox to simulate the effects of noise sources on decoherence
#  
# Author : Christophe Valahu 
# Date : 28/04/20
#
#######################################################################


import numpy as np
from math import pi, sqrt

def opamp_spectrum(en_white, en_flicker) :
    ''' Characteristic op-amp noise spectrum 
        
    Parameters
    ----------
    en_white : float   
        white noise level of PSD
    
    en_flicker : float
        1/f noise level of PSD
        
    
    '''
    return 1/en_flicker + en_white

def macro_voltage_spectrum(w) :
    ''' Example PSD for voltage noise on the mascroscopic experiment. 
    
    Parameters
    ----------
    w : float  
        angular frequency
    
    '''
    def Sj(w) :
    #Thermal Noise
        G_0 = 8.19e-9
        W_C = 2*pi*6.5e3
        N = .34
        K_B = 1.38064852e-23
        return 200*((G_0*(1+(6.28*w/W_C)**(2*N))**(-1.5))**2) + 4*K_B*300*0.4
    
    def Ss(w) :
    # Noise out of Voltage generator (HV 40-16)
        EN_FLICKER = (2.5e-6)**2
        EN_WHITE = (3.3e-9)**2
        if w < 0.1 :
            return EN_FLICKER*2*pi/0.1 + EN_WHITE
        else :
            return EN_FLICKER * 2 * pi/w + EN_WHITE
        
    def Filter(w, wc, n) :
    #Low pass butterworth filter
        return 1 / (1 + (w/wc)**(2*n))
    
    def Gauss(w, w0, sigma, A) :
    # Noise source modeled as a gaussian
        return (A * np.exp(-(w - w0)**2/(2*sigma**2)))**2
        
    # Filter parameters
    FILTER_WC = 2*pi*30.5
    FILTER_ORDER = 2
    
    # Gaussian noise source parameters
    GAUSS_W0 = 2*pi*20e3
    GAUSS_A = 7e-9
    GAUSS_SIGMA = 2*pi*4e3
    
    return abs(Ss(w)*Filter(w, FILTER_WC, FILTER_ORDER) + Sj(w) + Gauss(w, GAUSS_W0, GAUSS_SIGMA, GAUSS_A))