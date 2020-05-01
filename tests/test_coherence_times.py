import sys, os
sys.path.insert(0, os.path.abspath('../src/'))
import quantumnoise as qn
import pytest
import numpy as np

    
_test_data_source = [
    qn.Noise([qn.Source(lambda w: 10e-24, 1, 0)]),
    qn.Noise([qn.Source(lambda w: 10e-24, 0, 1)]),
    qn.Noise([qn.Source(lambda w: 10e-24, 1, 1)]),
    qn.Noise([qn.Source(lambda w: (10e-12)/w, 1, 1)])
]

_test_data_t1 = [ 
    2.05796185217845 , 
    np.inf, 
    2.05796185217845, 
    1.2930555672343727e-7
]

@pytest.mark.parametrize("noise, t1_expected", zip(_test_data_source, _test_data_t1))
def test_T1(noise, t1_expected) :
    
    Omw = 2*np.pi*10e3
    
    t1 = noise.T1(Omw)
    
    assert np.allclose(t1, t1_expected, atol = 1e-6)

    
_test_data_tphi = [ 
    2.05796185217845 , 
    np.inf, 
    2.05796185217845, 
    1.2930555672343727e-7
]

_test_data_N = [ 
    0 , 
    0, 
    1, 
    6
]
    
@pytest.mark.parametrize("noise, t1_expected", zip(_test_data_source, _test_data_phi, _test_data_N))
def test_Tphi(noise, tphi_expected, N) :
    
    Omw = 2*np.pi*10e3
    
    tphi = noise.Tphi(Omw, N)
    
    assert np.allclose(tphi, tphi_expected, atol = 1e-6)
    
    
    
    
    
    