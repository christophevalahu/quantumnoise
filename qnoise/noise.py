# QUANTUMNOISE
# Toolbox to simulate the effects of noise sources on decoherence
#  
# Author : Christophe Valahu 
# Date : 28/04/20
#
#######################################################################

import matplotlib.pyplot as plt
import time
import numpy as np
from math import sqrt, pi, cos
import cmath
from scipy import interpolate, integrate, optimize
import warnings

class Source(object) :
    ''' Class for a noise source
    
    Parameters
    ----------
    noise_func : function
        function of the noise PSD, its only input is the angular frequency w
        
    Dx : float
        sigma_x coupling strength of the qubit to its noise bath. 
    
    Dz : float
        sigma_z coupling strength of the qubit to its noise bath.
        
    sigma : float
        Standard deviation of the noise PSD
        
    Methods 
    --------
    show(w_min, w_max)
        Plots the noise spectrum
        
    '''
    def __init__(self, noise_func = None, Dx = 1, Dz = 1, sigma = 0) :
        self.noise_func = noise_func
        self.Dx = Dx
        self.Dz = Dz
        
        try :
            self.sigma = self.__sigma__()
        except :
            self.sigma = np.inf
            warnings.warn('Warning : Integral of spectrum diverges!')
            
        return None
    
    def __sigma__(self) :
        return sqrt(2*integrate.quad(lambda w: self.noise_func(w), 0.01, np.inf)[0])
    
    def show(self, w_min = 2*pi*0.1, w_max = 2*pi*1e5) : 

        fig = plt.figure()
        ax = fig.add_subplot(111)
        w_list = np.linspace(np.log10(w_min/(2*pi)), np.log10(w_max/(2*pi)), 100)
        PSD_list = [self.noise_func(2*pi*(10**w)) for w in w_list]
        w_list = [(10**w) for w in w_list]
        plt.loglog(w_list, PSD_list)
        plt.ylim(min(PSD_list)*0.5, max(PSD_list)*1.5)
        plt.xlim(w_min/(2*pi), w_max/(2*pi))
        plt.grid()
        plt.xlabel(r'$\omega/2\pi \ [Hz]$')
        plt.ylabel(r'$PSD$')
        plt.show()

        
        
class Noise(object) :
    ''' Class for decoherence due to noise
    
    Parameters
    ----------
    sources : list of Source
        A list of noise sources which contribute to decoherence. 
    
    Methods 
    --------        
    fdecay(t, Omw, N)
        Decay function for a free evolving qubit. Omw is the dressed state splitting
        and N is the number of refocussing pi pulses.
        
    fdecayrabi(t, OmR, OmW)
        Decay function for a driven evolution. OmR is the rabi frequency,
        Omw is the dressed state splitting.
        
    Tphi(Omw, N)
        Dephasing time
    
    T2(Omw, N)
        Coherence time
        
    T1(Omw)
        Lifetime
        
    '''
    
    def __init__(self, sources = None):
        
        if sources == None :
            self.sources = []
        else :
            if not isinstance(sources, list) or not all([isinstance(s, Source) for s in sources]) :
                raise TypeError('Input must be a list of Source objects.')
            else :
                self.sources=  sources
        
    
    def __gfilter__(self, w, t, N = 0, tpi = 0 ) :
        return 1/((w * t)**2) * abs(1 + 
            ((-1)**(N+1))*cmath.exp(1j * w * t)+
            2*sum([((-1)**(i+1))*cmath.exp(1j*w*t/(N+1)*(i+1))*cos(w*tpi/2) for i in range(N)]))**2
    
    def __fz__(self, t, N) :
        return  np.exp(- (t**2)* \
                    sum([source.Dz**2 *integrate.quad(lambda w: source.noise_func(w)*self.__gfilter__(w, t, N), 
                                                  0.01, np.inf, maxp1 = 200)[0] for source in self.sources]))

    def show(self, w_min = 2*pi*0.1, w_max = 2*pi*1e6) : 

        fig = plt.figure()
        ax = fig.add_subplot(111)
        w_log_list = np.linspace(np.log10(w_min/(2*pi)), np.log10(w_max/(2*pi)), 200)
        w_list = [(10**w) for w in w_log_list]
        ymax = 0
        ymin = 0
        
        total_PSD = np.zeros(200)
        
        for source in self.sources:
            PSD_list = [source.Dx**2 * source.noise_func(2*pi*(10**w)) for w in w_log_list]
            plt.loglog(w_list, PSD_list)
            ymax = max([max(PSD_list),ymax])
            ymin = min([min(PSD_list),ymax])
            
            total_PSD = [sum(x) for x in zip(total_PSD, PSD_list)]
        
        plt.scatter(w_list, total_PSD, s=5, facecolors='none', edgecolors='k')
        
        plt.ylim(ymin*0.5, ymax*1.5)
        plt.xlim(w_min/(2*pi), w_max/(2*pi))
        plt.grid()
        plt.xlabel(r'$\omega/2\pi \ [Hz]$')
        plt.ylabel(r'$S_B \  [T^2/Hz]$')
        plt.show()
    
    def fdecay(self, t, Omw, N = 0) :
        
        return self.__fz__(t, N) * np.exp(- t/(2*self.T1(Omw)))
    
    def fdecayrabi(self, t, OmR, Omw) :
        
        def _spectrum_dz(w) :
            return sum([source.Dz**2 * source.noise_func(w) for source in self.sources])
        
        sigma_total = sqrt(2*integrate.quad(_spectrum_dz, 0.01, np.inf)[0])
       
        def _chi(t) :
            return (1+((sigma_total**2)/OmR*t)**2)**(-1/4)
        
        Gammav = _spectrum_dz(OmR)/2
        
        GammaR = 3/4 /self.T1(Omw) + 1/2 * Gammav
        
        return _chi(t) * np.exp(-t*GammaR)
    
    def Tphi(self, Omw, N = 0) :
    
        def __minfun__(t) :
            return abs(np.exp(-1) - self.__fz__(10**t, N))
        
        bnds = [[-6, 2]]
        Tphi_0 = -3
        res = optimize.minimize(__minfun__, Tphi_0, method='Nelder-Mead', tol=1e-10, bounds = bnds)
        
        if abs(res.fun) >= 1e-6 :
            warnings.warn('A solution was not found : Scipy.optimize.minimize did not converge.')
            
        return 10**res.x[0]
        
    def T2(self, Omw, N = 0) :
        
        def __minfun__(t) :
            return abs(np.exp(-1) - self.fdecay(10**t, Omw, N)) 
        
        bnds = [[-6, 2]]
        T2_0 = -3
        res = optimize.minimize(__minfun__, T2_0, method='Nelder-Mead', tol=1e-10, bounds = bnds)
        
        if abs(res.fun) >= 1e-6 :
            warnings.warn('A solution was not found : Scipy.optimize.minimize did not converge.')
        
        return 10**res.x[0]
    
    def T1(self, Omw) :
        MU_B =  9.274009994e-24
        HBAR = 1.0545718e-34

        sb_tot = sum([source.Dx**2 * source.noise_func(Omw) for source in self.sources])
        
        if sb_tot == 0 :
            return np.inf
        else :
            return 1/(2*pi)*(HBAR/MU_B)**2 /sb_tot
       

        