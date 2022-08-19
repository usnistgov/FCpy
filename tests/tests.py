# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 13:40:17 2022

@author: Evan Groopman
@affiliation: National Institute of Standards and Technology (NIST)
@email: evan.groopman@nist.gov
"""

import pytest
import numpy as np

from FCpy.FCpy import FCpy as FC

def test_FC_poisson():
    assert np.all(FC.FC_poisson(0,0,0, conf=0.95, useCorrection= False).round(3) == np.array([0.   , 3.092]))
    
def test_FC_gauss():
    assert np.all(FC.FC_gauss(1, conf=0.95).round(3) == np.array([0.   , 2.969]))
    
def test_FC_poisson123():
    assert np.all(np.array(FC.FC_poisson123(0,0,0,0.95)).flatten().round(3) == np.array([0.   , 1.292, 0.   , 3.092, 0.   , 5.886]))

# if __name__ == '__main__':
#     pytest.main()