##########################
Physics of Decoherence
##########################

There are two sources of decoherence: depolarization and dephasing. 
They each have an associated decay time: lifetime and dephasing time. Together, both decays
cause decoherence. 

**************
Depolarization
**************

The depolarization rate (inverse of the lifetime) is defined as :

.. math:: 
	\Gamma_1 = \frac{2 \pi \mu_B^2}{\hbar^2} S_B(\omega)

	
where '$S_B(\omega)$' is the Power Spectral Density of magnetic field noise. :math:'\omega' corresponds to the frequency splitting
of the qubit. For dressed states we consider leakage to spectator state, which is analogous to depolarization, and so :math:'\omega = \frac{\Omega_{\mu w}}{\sqrt{2}}'.

*********
Dephasing
*********

The effects of dephasing strongly depend on the shape of the PSD. In the simplest case, if the PSD is constant accross the spectrum, we find that :

.. math::

	\Gamma_\phi = \pi S_{\delta\omega_z}(\omega = 0)
	:label: eq:Gammaphi
	

For an arbitrary noise spectrum, it is more difficult to derive an analytic expression for the dephasing rate. One can instead derive decay functions, which describe
the decay of fidelity over time. The decay function of dephasing is :

.. math::

	f_{\phi}(t) = exp(-\frac{t^2}{2}\int^{+\inf}_{-\inf} S_{\delta\omega_z}g_N(\omega, t) d\omega)
	:label: eq:fphi
	
where :math:'g_N(t)' is a filter function. During the evolution of a qubit, one can apply N refocussing pulses,
which decouple the qubit from low frequency noise. For N>0, the filter function acts as a band pass filter who's
center frequency increases with N. Therefore, for a 1/f spectrum, it is advantageous to increase the number of 
refocussing pulses. The filter function :math:'g_N(t)' is defined as :


.. math::
	
	g_N(\omega, t) = \frac{1}{(\omega t)^2}|1 + (-1)^N e^{i \omega t} + 2 \sum^N_{j=1}(-1)^j e^{i\omega j \frac{t}{N+1}}cos(\omega T_\pi/2)|^2
	:label: eq:gN
	
where Tpi is the time of a single refocussing pulse.

In order to obtain a characteristic dephasing time, one can solve the following equation:

.. math::
	t^2 = 2 \int^{+\inf}_{-\inf} S_{\delta\omega_z}g_N(\omega, t) d\omega

***********	
Decoherence
***********

Decoherence is the simultaneous effects of depolarization and dephasing. For a constant PSD, the decoherence rate is defined as:

.. math::
	
	\Gamma_2 = \frac{1}{2}\Gamma_1 + \Gamma_\phi
	:label: eq:Gamma2

However, as is the case for dephasing, most PSD contain colored noise. In such a case, the decay due to decoherence is a linear
combination of dephasing and depolarization decays. The decay function is :

.. math::
	f_2(t) = f_\phi(t)exp(-\Gamma_1 t/2) = exp(-\frac{t^2}{2}\int^{+\inf}_{-\inf} S_{\delta\omega_z}g_N(\omega, t) d\omega - \Gamma_1 \frac{t}{2})
	:label: eq:f2

Once again, the decay is not purely exponential, it is therefore difficult to obtain a analaytical expression for the decoherence time. However,
a characteristic decoherence time can be extracted by solving :math:'f_2(t) = e^{-1}'
	
	
****************
Driven evolution
****************

A final scenario arises when a weak field is resonant with the qubit transition frequency. This would be the case of a 
single qubit gate, which rotates with a frequency equal to the Rabi frequency :math:'\Omega_R'. There exists a first 
decay due to quasi-static noise, which is given by :

.. math::
	
	\chi(t) = (1 + (\sigma_{\delta\omega_z}t)^2)^{-1/4},
	
	\sigma_{\delta\omega_z}^2 = \int^{+\inf}_{-\inf} S_{\delta\omega_z}d\omega 
	:label: eq:chi
	
	
A second decay arises from noise resonant with the Rabi frequency. This decay is given by :

.. math::
	
	\Gamma_R = \frac{3}{4}\Gamma_1 + \frac{1}{2}\Gamma_\nu,
	:label: eq:GammaR
	
The final decay function can therefore be found to be :

.. math::
	
	f_R(t) = \chi(t) e^{-\Gamma_R t}
	:label: eq:fR
	
	
	
	
	