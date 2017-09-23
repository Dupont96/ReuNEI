from __future__ import print_function
import time
import pyatomdb
import pylab
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sunpy.map
import sunpy.data.sample
import pandas as pd


from pyatomdb.atomic import Ztoelsymb, Ztoelname
from pyatomdb.apec import solve_ionbal_eigen
from pandas import DataFrame

Aelements = ['H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
            'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
            'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U']

Belements = ['H', 'He', 'C', 
            'N', 'O', 'Ne',
            'Mg', 'Si', 'S', 
            'Ar', 'Ca', 'Fe']

data = pyatomdb.spectrum.Session()
#Ion Balance Calculations 


#Hydrogen-Helium Ratio
He_per_he = 0.1

for idx, val in enumerate(Belements):
    b = Ztoelname(idx+1)
    GenChargeStates = solve_ionbal_eigen(idx, 1e6)
    a = pd.DataFrame(data = GenChargeStates, columns=[b])
    print(a)
    a.to_csv('test.csv', index=True,header=True)

