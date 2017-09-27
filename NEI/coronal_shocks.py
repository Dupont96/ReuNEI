from __future__ import print_function

import numpy as np 
import numba as nb 
import matplotlib.pyplot as plt 
import drms 
import sunpy 
import sunpy.map 


c = drms.Client()


def something ():
    intensites = np.genfromtxt('aia_ascii.txt', dtype=None, delimiter=',')

    print(intensites)
something()