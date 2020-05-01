# QNOISE
# Toolbox to simulate the effects of noise sources on decoherence
#  
# Author : Christophe Valahu 
# Date : 28/04/20
#
#######################################################################

import numpy as np
from math import pi, sqrt
import warnings

        
def dw_db(state) :

    MU_B =  9.274009994e-24
    HBAR = 1.0545718e-34
    G_J = 2
    W_HF = 2*pi*12.642812118500e9
    B = 10e-4
        
    # dwdB is the change in frequency per change in magnetic field
    if state is "dressed" :
        return B/W_HF*((G_J*MU_B/HBAR)**2)/sqrt(1+(B*G_J*MU_B/(W_HF*HBAR))**2)/sqrt(2)
    elif state is "clock" :
        return B/W_HF*((G_J*MU_B/HBAR)**2)/sqrt(1+(B*G_J*MU_B/(W_HF*HBAR))**2)
    elif state is "bare":
        return 1/2*W_HF*(G_J * MU_B/(W_HF*HBAR) - B*G_J**2 * MU_B**2/((W_HF*HBAR)**2*sqrt(1+(B*G_J*MU_B/W_HF/HBAR)**2)))
    else :
        warnings.warn('State was not recognized!')
        return 0
                

def voltage_Dx(dzB, v, d, eta, Omw = 0) :
    ''' Sigmax coupling of voltage noise
    
    Parameters
    ----------
    dzB : float
        Magnetic field gradient [T/m]
    
    v : float
        secular frequency [Hz]
    
    d : float
        distance to electrodes [m]
    
    eta : float 
        geometric factor [a.u.]
    
    Omw : float 
        dressed state splitting [Hz]
    
    Notes
    ---------
    Omw can be set to 0 when v >> Omw, which is usually the case
    
    '''
    m =  171*1.6726219e-27;
    q = 1.60217662e-19;
    return q*dzB*eta/(m*(v**2 - Omw**2)*d)
      
    
    
def voltage_Dz(dzB, v, d, eta, Omw = 0, state = "dressed") :
    ''' Sigmaz coupling of voltage noise
    
    Parameters
    ----------
    dzB : float
        Magnetic field gradient [T/m]
    
    v : float
        secular frequency [Hz]
    
    d : float
        distance to electrodes [m]
    
    eta : float 
        geometric factor [a.u.]
    
    Omw : float 
        dressed state splitting [Hz]
    
    state : string, ["dressed", "clock", "bare"]
        the qubit's spin state 
    
    Notes
    ---------
    Omw can be set to 0 when v >> Omw, which is usually the case
    
    '''
    
    Dx = voltage_Dx(dzB, v, d, eta, Omw)
    
    # Dx can be written as dB/dV, ie change in magnetic field per change in voltage
    # noise. We can therefore use dw/dV = dw/dB * dB/dV = dw/dB * Dx
    return dw_db(state) * Dx
    

def current_Dx(dx, z0) :
    ''' Sigmax coupling of current noise
    
    Parameters
    ----------
    dx : float
        Distance of the ion from the magnetic gradient nil [m]
    
    z0 : float
        distance of the ion from the chip [m]
    '''
    MU_0 =  1.25663706212e-6 #T m/A
    
    return dx * 9 * MU_0/(8*pi * sqrt(3) * (z0**2))
    
def current_Dz(dx, z0, state = "dressed") :
    ''' Sigmaz coupling of current noise
    
    Parameters
    ----------
    dx : float
        Distance of the ion from the magnetic gradient nil [m]
    
    z0 : float
        distance of the ion from the chip [m]
        
    state : string, ["dressed", "clock", "bare"]
        the qubit's spin state 
    '''
    
    return dw_db(state) * current_Dx(dx, z0)
    
