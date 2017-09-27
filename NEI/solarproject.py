from __future__ import print_function

import time
import pyatomdb
import pylab
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sunnei
import glob
import xarray as xr
import sys
import ChiantiPy.core as ch
import os.path
import operator
import matplotlib
import pprint

from astropy.io import fits
from pyatomdb.atomdb import get_data, lorentz_neicsd, lorentz_power
from pyatomdb.atomic import Ztoelsymb, Ztoelname
from pyatomdb.apec import *
from pandas import DataFrame
from sunnei import read_atomic_data
from astropy.io import fits
from numpy import linalg as LA
from scipy import interpolate

sys.path.insert(0, '/home/mdupont/2017_Projects/SunNEI/sunnei/applications/')

sys.path.insert(0, '/home/mdupont/2017_Projects/ionica/ReuNEI/sunnei/applications/')

import cmeheat2
import cmeheat


Zlist = pd.Series(np.arange(28)+1,
            index=['H','He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
            'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni'])


Belements = ['H', 'He', 'C', 
            'N', 'O', 'Ne',
            'Mg', 'Si', 'S', 
            'Ar', 'Ca', 'Fe']

data = pyatomdb.spectrum.Session()

pp = pprint.PrettyPrinter(indent=3)
def dummy_plot():
    logT = np.linspace(4.0,8.0,161) # Delta T = 0.025
    
    T = 10**logT

    Lambda = np.array(
        [0.662E+00,0.113E+01,0.190E+01,0.312E+01,0.492E+01,0.725E+01,
         0.968E+01,0.115E+02,0.123E+02,0.127E+02,0.123E+02,0.116E+02,
         0.107E+02,0.988E+01,0.932E+01,0.904E+01,0.901E+01,0.934E+01,
         0.991E+01,0.107E+02,0.115E+02,0.125E+02,0.137E+02,0.151E+02,
         0.168E+02,0.188E+02,0.213E+02,0.241E+02,0.275E+02,0.313E+02,
         0.357E+02,0.406E+02,0.460E+02,0.518E+02,0.576E+02,0.628E+02,
         0.671E+02,0.699E+02,0.710E+02,0.702E+02,0.677E+02,0.633E+02,
         0.580E+02,0.526E+02,0.480E+02,0.445E+02,0.422E+02,0.407E+02,
         0.397E+02,0.390E+02,0.382E+02,0.375E+02,0.368E+02,0.363E+02,
         0.357E+02,0.348E+02,0.333E+02,0.311E+02,0.280E+02,0.246E+02,
         0.214E+02,0.187E+02,0.166E+02,0.152E+02,0.142E+02,0.136E+02,
         0.133E+02,0.132E+02,0.130E+02,0.128E+02,0.124E+02,0.120E+02,
         0.116E+02,0.112E+02,0.110E+02,0.109E+02,0.109E+02,0.111E+02,
         0.112E+02,0.113E+02,0.114E+02,0.115E+02,0.115E+02,0.114E+02,
         0.113E+02,0.112E+02,0.111E+02,0.109E+02,0.106E+02,0.102E+02,
         0.974E+01,0.914E+01,0.840E+01,0.763E+01,0.687E+01,0.616E+01,
         0.553E+01,0.500E+01,0.457E+01,0.423E+01,0.396E+01,0.376E+01,
         0.360E+01,0.350E+01,0.344E+01,0.341E+01,0.341E+01,0.345E+01,
         0.350E+01,0.357E+01,0.364E+01,0.371E+01,0.377E+01,0.381E+01,
         0.382E+01,0.379E+01,0.371E+01,0.359E+01,0.343E+01,0.325E+01,
         0.305E+01,0.286E+01,0.267E+01,0.249E+01,0.234E+01,0.221E+01,
         0.210E+01,0.201E+01,0.193E+01,0.187E+01,0.183E+01,0.179E+01,
         0.177E+01,0.175E+01,0.174E+01,0.173E+01,0.173E+01,0.173E+01,
         0.174E+01,0.176E+01,0.177E+01,0.179E+01,0.182E+01,0.184E+01,
         0.187E+01,0.191E+01,0.194E+01,0.198E+01,0.202E+01,0.206E+01,
         0.210E+01,0.214E+01,0.219E+01,0.224E+01,0.229E+01,0.234E+01,
         0.239E+01,0.244E+01,0.250E+01,0.256E+01,0.261E+01,]
        )*1e-23
    
    assert T.size == Lambda.size, 'Mismatch in cooling function tables'

    lol = interpolate.interp1d(T, Lambda, fill_value=0.0)
    return lol

#P_cool coefficients
def lorentz_power(version):
  """
  Calculate the power emitted from 13.6eV to 13.6keV in a 1m^3 slab of
  plasma with n_e=1e6m^-3.

  Parameters
  ----------
  version : string
    The version string

  Returns
  -------
  None

  """
  
  # get the elemental abundances
  ag89=pyatomdb.atomdb.get_abundance(abundset='AG89')

  # open results file
  f = open('Power_atomdb_%s.dat'%(version),'w')

  # specify energy bins from 0.01 to 100 keV
  ebins = numpy.linspace(0.01, 100.0, 100001)
  
  # calculate the central energy of each bin in ERG
  energy = pyatomdb.const.ERG_KEV*(ebins[1:]+ebins[:-1])/2
  s = "Z  z1 "


  # print out temperature vector
  for iT in range(51):
    flt = 4.0+(0.1*iT)
    s += "        %.1f"%(flt)

  f.write("%s\n"%(s))

  # open line and continuum emissivity files
  ldat = pyatomdb.pyfits.open('/home/mdupont/atomdb/apec_v3.0.8_nei_line.fits')
  cdat = pyatomdb.pyfits.open('/home/mdupont/atomdb/apec_v3.0.8_nei_comp.fits')


  # cycle through each ion
  for element in Belements:
    Z = Zlist[element]
    for z1 in range(Z+1):
      # holder for total energy at each temperature
      tot_e = numpy.zeros(51)
      
      
      s = "%2i %2i"%(Z, z1)
      
      for iT in range(2,53):
        # calculate spectrum in photons cm^3 s^-1 bin^-1, with AG89 abundance built in
        spec= pyatomdb.spectrum.make_ion_spectrum(ebins, iT, Z, z1, linefile = ldat,\
                                            cocofile = cdat,\
                                            nei=True)
      
        # convert to erg cm^3 s^-1 bin^-1
        e = spec*energy
        
        # sum the energy over the bins, then remove abundance
        tot_e[iT-2] = sum(e)/ag89[Z]
        
        # ensure minimum value of 1e-60 erg cm^s^-1
        
        s+= " %10.6f"%(numpy.log10(max(1e-60, tot_e[iT-2])))
        
      # print results for ion
      f.write("%s\n"%(s))
  f.close()

def get_ion_fractions():
    f = open('ion_fractions.txt', 'w')
    nei_frac = cmeheat2.cmeheat_track_plasma2()
    l = "Z z1"

    f.write('%s\n'%(l))

    for element in Belements:
        Z = Zlist[element]
        for z1 in range(Z+1):
            l = "%2i %2i"%(Z,z1)
            ion_frac = nei_frac['ChargeStates'][element][-1][z1]
            l+="  %g"%(ion_frac)

            f.write("%s\n"%(l))
    f.close()

def ionica():
    nei_frac = np.loadtxt('ion_fractions.txt', skiprows=1)
    pp = pprint.PrettyPrinter(indent=3)

    
    frac_dict = {}
    for i, j in enumerate(nei_frac):
        zed, ion = nei_frac[i][:2]
        frac = nei_frac[i][2:]

        frac_dict[zed, ion] = {
            'ion_fraction':frac
        }
    ion_frac = {}
    i=0
    
    
    for element in Belements:
        Z = Zlist[element]
        ncharge = Z+1
        
        ion_frac[Z] = {}

        for z1 in range(ncharge):
            ion_frac[Z][z1] = frac_dict[Z, z1]['ion_fraction']
    
    
    return ion_frac


def gen_cooling():
    cool = os.path.exists('/home/mdupont/2017_Projects/ionica/NEI/Power_atomdb_3.0.8.dat')
    nei_calc  = os.path.exists('ion_fractions.txt')
    matplotlib.rcParams.update({'font.size': 25})
    f = open('sigma2.txt', 'w')
    if not cool:
        lorentz_power('3.0.8')
        p_cool = np.genfromtxt('Power_atomdb_3.0.8.dat', dtype=None, delimiter=',') #our cooling coefficient (Power)
    else:
        p_cool = np.genfromtxt('Power_atomdb_3.0.8.dat', dtype=None , delimiter=',')
    abund = pyatomdb.atomdb.get_abundance(abundset='AG89')

    log_matrix = np.zeros((153,53), dtype=float) #Matrix used to store our values for the cooling rates

    for idx, line in enumerate(p_cool):
        if line[idx] != 'Z':
            log_matrix[idx] = line.split()[:]
        else:
            Te = np.array(line.split()[2:], dtype=float)
            
    LogT = np.linspace(Te[0], Te[-1], 152)
    T = 10**LogT
    Te = 10**Te
    c_matrix = 10**log_matrix[1:] #Scales the cooling terms which were initially in Log base 10
    Te_vec = Te.reshape(1,-1)


    rate_matrix = log_matrix[1:] #truncate matrix to exclude extraneous array of 0s

    Lambda ={} #intialize the atoms dictionary
    for idx, val in enumerate(rate_matrix):
        zed, ion = rate_matrix[idx][:2] #The corresponding ion for each Z
        cool_terms = c_matrix[idx][2:] #The cooling_terms array that thus corresponds to each Z and z1
        for k, temp in enumerate(Te):

            Lambda[zed, ion, temp] = {
                'rate': cool_terms[k]
            }
    
    
    pp.pprint(Lambda)
    

    
    '''
    sum_list = []
    r= []
    nei_frac = ionica()

    s = "Z  z1"

    for iT in range(51):
       flt = 4.0+(0.1*iT)
       s += "        %.1f"%(flt)

    f.write("%s\n"%(s))
    
    
    for el in Belements:
        Z = Zlist[el]
        ncharge = Z+1
        for z1 in range(ncharge):
            s = "%2i %2i"%(Z, z1)
            lul  = [x*nei_frac[Z][z1]*abund[Z] for x in b]
            for item in lul:
                
                s += "  %10.3g"%(item)
            f.write(" %s\n"%(s))
            sum_list.append(lul)
            
            
    f.close()
    arr = np.array(sum_list)
    plip = np.sum(sum_list)
    np.savetxt('nothing2.txt', arr)
    '''
    
def get_lambda(Te,Z,z1):
    cool = os.path.exists('Power_atomdb_3.0.8.dat')
    if not cool:
        lorentz_power('3.0.8')
        p_cool = np.genfromtxt('Power_atomdb_3.0.8.dat', dtype=None, delimiter=',')
    else:
        p_cool = np.genfromtxt('Power_atomdb_3.0.8.dat', dtype=None, delimiter=',')
    log_matrix = np.zeros((153,53), dtype=float)

    
    for i, line in enumerate(p_cool):
        if line[i] != 'Z':
            log_matrix[i] = line.split()[:]
        else:
            LogT = np.array(line.split()[2:], dtype=float)
    log_rate = log_matrix[1:]

    temp = 10**LogT
    c_matrix = 10**log_matrix[1:]

    #Declare Lambda dictionary storing the respective cooling rates per element-ion pair
    Lambda = {}
    gamm = {}
    for idx, val in enumerate(log_rate):
        zed, ion = log_rate[idx][:2]
        cooling_terms = c_matrix[idx][2:]

        for k, T in enumerate(temp):

            Lambda[zed, ion, T] = {
                'rate': cooling_terms[k]
            }
    return Lambda[Z,z1, Te]['rate']



def get_cool(zed, ion):
    cool = os.path.exists('Power_atomdb_3.0.8.dat')
    if not cool:
        lorentz_power('3.0.8')
        p_cool = np.genfromtxt('Power_atomdb_3.0.8.dat', dtype=None, delimiter=',')
    else:
        p_cool = np.genfromtxt('Power_atomdb_3.0.8.dat', dtype=None, delimiter=',')
    log_matrix = np.zeros((153,53), dtype=float)

    
    for i, line in enumerate(p_cool):
        if line[i] != 'Z':
            log_matrix[i] = line.split()[:]
        else:
            LogT = np.array(line.split()[2:], dtype=float)
    log_rate = log_matrix[1:]

    T = 10**LogT
    c_matrix = 10**log_matrix[1:]

    #Declare Lambda dictionary storing the respective cooling rates per element-ion pair
    gamma = {}
    for idx, val in enumerate(log_rate):
        Z, z1 = log_rate[idx][:2]
        cooling_terms = c_matrix[idx][2:]
        gamma[Z,z1] = {
            'rates': cooling_terms
        }

    return gamma[zed, ion]['rates']

    
def yes():
    matplotlib.rcParams.update({'font.size':25})
    i = 1
    Lambda = get_cool
    abund = pyatomdb.atomdb.get_abundance(abundset='AG89')
    a = []
    Te = np.logspace(4.0,9.0,51)
    nei_frac = ionica()

    logT = np.linspace(4.0,8.0,161) # Delta T = 0.025
    
    T = 10**logT

    for el in Belements:
        Z = Zlist[el]
        ncharge = Z+1
        for z1 in range(ncharge):
            calc = Lambda(Z,z1)*abund[Z]*nei_frac[Z][z1]
            a.append(calc)
    a = np.array(a)

    b = (a.sum(axis=0))/153
    
    F = interpolate.interp1d(Te, b, fill_value=0.0)

    gamma = dummy_plot()

    
   
    ax = plt.subplot(211)
    plt.loglog(Te, F(Te), '-x', label='New')
    plt.ylabel('$\Lambda(T, Z,z1)$')
    plt.title('New VS Old')

    plt.subplot(212, sharex=ax)
    plt.loglog(T, gamma(T), '-x', label='Old')
    plt.ylabel('$\Lambda (Equilibrium)$')
    plt.xlabel('Temperature in K')

    plt.legend()
    plt.show()








def ion_rec():
    matplotlib.rcParams.update({'font.size': 30})
    #f = open('recomb_rates_atomdb.txt','w')
    #f2 = open('recomb_rates_chianti.txt', 'w')

    

    
    d={}
    t=np.logspace(5,9, 100)
    ith =np.arange(2,11)
    ion_list = []
    rec_list = []
    for i in ith:
        ion, rec = pyatomdb.atomdb.get_ionrec_rate(t, False, Te_unit = 'K',Z=26,z1=i,datacache=d)
        ion_list.append(ion)
        rec_list.append(rec)

        #ax1 = plt.subplot(211)
        #plt.plot(t, rec, label='Fe %s +'%i, linestyle='-.')
        #plt.loglog(t, rec)
        #plt.ylabel('Recombination Rates')
        
        #ax2=plt.subplot(212)
        #plt.plot(t, ion, label='Fe %s +'%i, linestyle='-.')
        #plt.loglog(t, ion)
        #plt.ylabel('Ionization Rates')
    ioniz = np.array(ion_list)
    recomb = np.array(rec_list)

    fe_x = ['fe_3','fe_4','fe_5','fe_6','fe_7','fe_8','fe_9','fe_10','fe_11']
    c_list=[]
    irates = []
    rrates = []
    for num, ion in enumerate(fe_x):
        fex = ch.ion(ion,temperature =t)
        fex.ionizRate()
        fex.recombRate()
        plop = fex.IonizRate
        plip = fex.RecombRate
        c_list.append(fex)
        irates.append(plop['rate'])
        rrates.append(plip['rate'])

        #plt.subplot(211)
        #plt.loglog(t, plip['rate'], label='Fe %s +'%(num+2))

        #plt.subplot(212)
        #plt.loglog(t, plop['rate'], label='Fe %s +'%(num+2))
        #plt.xlabel('Temperature in K')

    irates = np.array(irates)
    rrates = np.array(rrates)

    print(rrates)
    print()
    print(recomb)

    for j, rate in enumerate(rrates):
        ax1 = plt.subplot(211)
        plt.loglog(t,rate,'-.')
        plt.loglog(t,recomb[j])
        plt.ylabel('Recomb Rates')
        plt.xlim(t[0],t[-1])

        ax2 = plt.subplot(212,sharex=ax1)
        plt.loglog(t, recomb[j]/rate, label='Fe %i +'%(j+2))
        plt.ylabel('Ratio AtomDB/CHIANTI')
        plt.xlabel('Temperature (K)')
        


   
    #plt.legend(loc='center left',bbox_to_anchor=(1,0.3), ncol=2)
    plt.legend(loc='upper right',ncol=2)
    plt.show()

ion_rec()



def mean_charge_state():
    pass

#For every element of charge Z, there are Z+1 possible ion states
def ion_states():
    for zed in Belements:
        nuc = Zlist[zed]
        Z_ion = nuc + 1
        name = Ztoelname(nuc)
        print(name,"yields",Z_ion,"possible charge states")


He_per_H = 0.1
#Hydrogen-Helium Ratio
He_per_H = 0.1


#Reading in AtomDB atomic data of high abundance elements
stuff=0
def ion_info():
    stuff={}
    for element in Belements:
        Z = Zlist[element] #assigns element index from Zlist to Z
        ElName = Ztoelname(Z)
        ElSymb = Ztoelsymb(Z)
        stuff[element] = {
            'Z':Z,
            'ElName':ElName,
            'ElSymb':ElSymb
        }
        for i in range(Z+1):
            if i != 0:
                print("The corresponding charge states of",ElName,"are",ElSymb,"+",i)
            else:
                print("Neutral",ElName)

def tabular_data():
    fits_list =[]
    for fit in glob.glob('/home/mdupont/atomdb/APED/*/*/*.fits'):
        f = fits.open(fit)
        tbdata = f[1].data
        fits_list.append(f)
        df = xr.DataArray(tbdata)
        print(df)

def data_grab():
    for element in Belements:
        Z = Zlist[element]
        for i in range(Z+1):
            z1 = i
            get_data(Z,z1,'IR')
            get_data(Z,z1,'DR')
        


def electron_density():
    zed = []
    for el in Belements:
        zed.append(Zlist[el])

    for idx, val in enumerate(zed):
        Ion = val, idx
    print(Ion, idx)

def compare():
    result_atomdb = cmeheat2.cmeheat_track_plasma2()
    result_chianti = cmeheat.cmeheat_track_plasma()
    elements = ['C', 'Si', 'Fe']
    for el in elements:
        nstates = Zlist[el] + 1
        for z1 in range(nstates):
            abund_atomdb = result_atomdb['ChargeStates'][el][-1][z1]
            abund_chianti = result_chianti['ChargeStates'][el][-1][z1]
            if not np.isclose(abund_atomdb,abund_chianti, rtol = 0.005, atol = 0.001):
                print(el, z1, abund_atomdb, abund_chianti)

class MyCache:







    def __init__(self):

        self.cache = {}
        self.max_cache_size = 1000




        def __contains__(self, key):
            return key in self.cache









    






        
                











