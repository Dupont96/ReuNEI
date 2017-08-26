from __future__ import print_function
import numpy as np
from scipy import interpolate
import os
import pyatomdb 
import matplotlib.pyplot as plt

def atomdb_cooling():
    cool = os.path.exists('/home/mdupont/2017_Projects/ionica/NEI/Power_atomdb_3.0.8.dat')
    if not cool:
        lorentz_power('3.0.8')
        p_cool = np.genfromtxt('../Power_atomdb_3.0.8.dat', dtype=None, delimiter=',') #our cooling coefficient (Power)
    else:
        p_cool = np.genfromtxt('../Power_atomdb_3.0.8.dat', dtype=None , delimiter=',')
    abund = pyatomdb.atomdb.get_abundance(abundset='AG89')
    chg_list = []
    som_matrix = np.zeros((153,53), dtype=float) #Matrix used to store our values for the cooling rates

    for idx, line in enumerate(p_cool):
        if line[idx] != 'Z':
            chg_list.append(np.array(line.split()[:], dtype=float))
            som_matrix[idx] = line.split()[:]
        else:
            Te = np.array(line.split()[2:], dtype=float)
            
    #nei_frac = cmeheat2.cmeheat_track_plasma2()
    
    rate_matrix = som_matrix[1:] #truncate matrix to exclude extraneous array of 0s

    Lambda ={} #intialize the atoms dictionary
    for idx, val in enumerate(rate_matrix):
        zed, ion = rate_matrix[idx][:2]
        cool_terms = rate_matrix[idx][2:] #The cooling_terms array that thus corresponds to each Z and z1
        
        Lambda [zed, ion] = {        #Creates a dictionary with tupled keys (Z,z1) which is called later
            'cooling_terms': cool_terms
        }
    sum_list = []
    for Z in range([1,2,6,7,8,9,10]):
        for z1 in range(Z+1):
            lol = sum(Lambda[Z,z1]['cooling_terms'])
            sum_list.append(lol)

    sum_l = np.array(sum_list)
    print(sum_l)

atomdb_cooling()  