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
from math import *
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
        self.sigma = self.__sigma__()
        return None
    
    def __sigma__(self) :
        return sqrt(2*integrate.quad(lambda w: self.noise_func(w), 0.01, np.inf)[0])
    
    def show(self, w_min = 2*pi*0.1, w_max = 2*pi*1e5) : 

        fig = plt.figure()
        ax = fig.add_subplot(111)
        w_list = np.linspace(log10(w_min/(2*pi)), log10(w_max/(2*pi)), 100)
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
        return  exp(- (t**2/4)* \
                    sum([source.Dz**2 *integrate.quad(lambda w: source.noise_func(w)*self.__gfilter__(w, t, N), 
                                                  0.01, np.inf, maxp1 = 200)[0] for source in self.sources]))
    def fz(self, t, N) :
        return  exp(- (t**2/4)* \
                    sum([source.Dz**2 *integrate.quad(lambda w: source.noise_func(w)*self.__gfilter__(w, t, N), 
                                                  0.01, np.inf, maxp1 = 200)[0] for source in self.sources]))
    
    def show(self, w_min = 2*pi*0.1, w_max = 2*pi*1e6) : 

        fig = plt.figure()
        ax = fig.add_subplot(111)
        w_log_list = np.linspace(log10(w_min/(2*pi)), log10(w_max/(2*pi)), 200)
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
        
        return self.__fz__(t, N) * exp(- t/(2*self.T1(Omw)))
    
    def fdecayrabi(self, t, OmR, Omw) :
        
        u = sum([(source.Dz**2)*(source.sigma**2) for source in self.sources])/OmR
        
        def chi(t) :
            return (1+(u*t)**2)**(-1/4)
        
        GammaR = 3/4 * 1/self.T1(Omw) + sum([1/2*(pi*source.noise_func(OmR)*source.Dz) for source in self.sources])
        
        return chi(t) * exp(-t*GammaR)
    
    def Tphi(self, Omw, N = 0) :
    
        def __minfun__(t) :
            return abs(exp(-1) - self.__fz__(10**t, N))
        
        bnds = [[-6, 2]]
        Tphi_0 = -3
        res = optimize.minimize(__minfun__, Tphi_0, method='Nelder-Mead', tol=1e-10, bounds = bnds)
        
        if abs(res.fun) >= 1e-6 :
            warnings.warn('A solution was not found : Scipy.optimize.minimize did not converge.')
            
        return 10**res.x[0]
        
    def T2(self, Omw, N = 0) :
        
        def __minfun__(t) :
            return abs(exp(-1) - self.fdecay(10**t, Omw, N)) 
        
        bnds = [[-6, 2]]
        T2_0 = -3
        res = optimize.minimize(__minfun__, T2_0, method='Nelder-Mead', tol=1e-10, bounds = bnds)
        
        if abs(res.fun) >= 1e-6 :
            warnings.warn('A solution was not found : Scipy.optimize.minimize did not converge.')
        
        return 10**res.x[0]
    
    def T1(self, Omw) :
        MU_B =  9.274009994e-24
        HBAR = 1.0545718e-34
      
        return 1/(2*pi)*(HBAR/MU_B)**2 /sum([source.Dx**2 * source.noise_func(Omw) for source in self.sources])
       

        
def voltage_Dx(dzB, v, d, eta, Omw = 0) :
    ''' Sigmax coupling of voltage noise
    
    Parameters
    ----------
    dzB : Magnetic field gradient
    
    v : secular frequency
    
    d : distance to electrodes
    
    eta : geometric factor
    
    Omw : dressed state splitting
    
    
    Notes
    ---------
    Omw can be neglected when v >> Omw, which is usually the case
    
    '''
    m =  171*1.6726219e-27;
    q = 1.60217662e-19;
    return q*dzB*eta/(m*(v**2 - Omw**2)*d)

def voltage_Dz(dzB, v, d, eta, Omw = 0, sensitive = False) :
    ''' Sigmax coupling of voltage noise
    
    Parameters
    ----------
    dzB : Magnetic field gradient
    
    v : secular frequency
    
    d : distance to electrodes
    
    eta : geometric factor
    
    Omw : dressed state splitting
    
    Sensitive : True if using magnetic sensitive states
    
    
    Notes
    ---------
    Omw can be neglected when v >> Omw, which is usually the case
    
    '''
    
    MU_B =  9.274009994e-24
    HBAR = 1.0545718e-34
    G_J = 2
    W_HF = 2*pi*12.642812118500e9
    B = 10e-4
        
    Dx = voltage_Dx(dzB, v, d, eta, Omw)
    
    # dwdB is the change in frequency per change in magnetic field
    if not sensitive :
        dwdB = B/W_HF*((G_J*MU_B/HBAR)**2)/sqrt(1+(B*G_J*MU_B/(W_HF*HBAR))**2)
    else:
        dwdB = 1/2*W_HF*(G_J * MU_B/(W_HF*HBAR) - B*G_J**2 * MU_B**2/((W_HF*HBAR)**2*sqrt(1+(B*G_J*MU_B/W_HF/HBAR)**2)))

    # Dx can be written as dB/dV, ie change in magnetic field per change in voltage
    # noise. We can therefore use dw/dV = dw/dB * dB/dV = dw/dB * Dx
    return dwdB * Dx

def v_noise_macro(w) :
    ''' Example PSD for voltage noise on the mascroscopic experiment. 

    Methods
    ------
    Sj : Thermal noise

    Ss : Noise out of the Voltage Generator (HV-40-16)

    Filter : Low pass butterworth filter

    Gauss : Noise source which was identified at 20 kHz, modeled by a gaussian
    
    '''
    def Sj(w) :
        G_0 = 8.19e-9
        W_C = 2*pi*6.5e3
        N = .34
        K_B = 1.38064852e-23
        return 200*((G_0*(1+(6.28*w/W_C)**(2*N))**(-1.5))**2) + 4*K_B*300*0.4
    
    def Ss(w) :
        EN_FLICKER = (3e-6)**2
        EN_WHITE = (3.3e-9)**2
        if w < 0.1 :
            return EN_FLICKER*2*pi/0.1 + EN_WHITE
        else :
            return EN_FLICKER * 2 * pi/w + EN_WHITE
        
    def Filter(w, wc, n) :
        return 1 / (1 + (w/wc)**(2*n))
    
    def Gauss(w, w0, sigma, A) :
        return (A * exp(-(w - w0)**2/(2*sigma**2)))**2
        
    # Filter parameters
    FILTER_WC = 2*pi*30.5
    FILTER_ORDER = 2
    
    # Gaussian noise source parameters
    GAUSS_W0 = 2*pi*20e3
    GAUSS_A = 7e-9
    GAUSS_SIGMA = 2*pi*4e3
    
    return abs(Ss(w)*Filter(w, FILTER_WC, FILTER_ORDER) + Sj(w) + Gauss(w, GAUSS_W0, GAUSS_SIGMA, GAUSS_A))