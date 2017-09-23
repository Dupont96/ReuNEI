from __future__ import print_function

import numpy as np 
import numba as nb 
import matplotlib.pyplot as plt 
import pandas as pd
import matplotlib
import time
import pyatomdb
import glob
import matplotlib
import scipy.integrate as integrate
import pprint
import datetime

from scipy import interpolate
from astropy.io import ascii, fits
from astropy.nddata import Cutout2D
from astropy import units as u
from scipy.optimize import curve_fit
from matplotlib import dates

position = (3500,1400)
size = (1000,1000)

pp = pprint.PrettyPrinter(indent=3)

matplotlib.rcParams.update({'font.size':25})

Zlist = pd.Series(np.arange(28)+1,
            index=['H','He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
            'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni'])


Belements = ['H', 'He', 'C', 
            'N', 'O', 'Ne',
            'Mg', 'Si', 'S', 
            'Ar', 'Ca', 'Fe']

#Constants 
h = 6.62607004e-27 #Plank's constant in ergs/s
cl = 2.998e10 #speed of light in cm/s
lam_171 = 171e-10 #wavelength in A
lam_193 = 193e-10 #wavelength in A
lam_211 = 211e-10 #wavelength in A
lam_335 = 335e-10 #wavelength in A
electrons = 1/(5.84794e-12) #ergs to electrons
gain_171 = 17.7 #electrons/DN 
gain_193 = 18.3 #electrons/DN
gain_211 = 18.3 #electrons/DN
gain_335 = 17.6 #electrons/DN

def spec_time():
    matplotlib.rcParams.update({'font.size': 25})
    arr =[]


    
    with open('aia_ascii.txt', 'r') as f:
        header1 = f.readline()
        header2 = f.readline()
        header3 = f.readline()
        header4 = f.readline()
        header5 = f.readline()

        columns = f.readline()

        for line in f:
            line = line.strip()
            vals = np.array(line.split(), dtype=float)
            arr.append(vals)
        
    
    
    arr = np.array(arr[:-1], dtype=float) #for some reason the file was adding an extra empty array at the very end of the list, so I cut it off
    time = arr[:, 0] #time[sec]
    band_171 = arr[:,1] #171 Angstrom band pass
    band_193 = arr[:,2] #192 Angstrom band pass
    band_211 = arr[:,3]*5 #211 Angstrom band pass scaled by 5
    band_335 = arr[:,4]*250 #335 Angstrom band pass scaled by 250

    
    f1 = interpolate.interp1d(time, band_171, fill_value=0.0)
    f2 = interpolate.interp1d(time, band_193, fill_value=0.0)
    f3 = interpolate.interp1d(time, band_211, fill_value=0.0)
    f4 = interpolate.interp1d(time, band_335, fill_value=0.0)

    hour = 5
    time_fmt = [] #Time format list
    for sec in time:
        m,s = divmod(1800+sec,60)
        h,m = divmod(m,60)
        time_fmt.append("%d:%02d:%02d"%(hour,m,s))
    tls = []

    for i, t in enumerate(time_fmt):
        dts = datetime.datetime.strptime(time_fmt[i],'%H:%M:%S')
        tls.append(dts)

    fds = dates.date2num(tls) #Converted times to numbers readable by matplotlib

    # matplotlib time format object
    hfmt = dates.DateFormatter('%H:%M:%S')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    major_ticks = np.arange(fds[0],fds[-1])

    ax.plot(fds, f1(time), label='171 Angstrom')
    ax.plot(fds, f2(time), label='193 Angstrom')
    ax.plot(fds, f3(time), label='211 Angstrom')
    ax.plot(fds, f4(time), label='335 Angstrom')
    

    ax.xaxis.set_major_locator(dates.MinuteLocator(interval=5))
    ax.xaxis.set_major_formatter(hfmt)
    ax.set_xlim(fds[0],fds[-1])

    #for label in ax.xaxis.get_ticklabels()[::2]:
     #   label.set_visible(False)

    plt.axvline(x=fds[45o   ],linestyle='--')
    plt.ylabel('Intensity DN s^-1 pix^-1')
    plt.xticks(rotation='vertical')
    plt.legend()
    plt.show() #Replicates Figure 2 from Ma et al. paper

    #return f1,f2,f3,f4, time
    return time_fmt, time, f4




def flux():
    wave_171 = {}
    wave_193 = {}
    wave_211 = {}
    wave_335 = {}

    #This block grabs the first set of sits files from roughly 5:39 UT (Time of shock appearance) til about 5:39:50
    b171_list = [fits.open('/home/mdupont/AIA_data/171_A/aia.lev1_euv_12s.20100613T0539'+n+'Z.171.image_lev1.fits',cache=True) for n in ['02','14','26','38','50']]
    b193_list = [fits.open('/home/mdupont/AIA_data/193_A/aia.lev1_euv_12s.20100613T0539'+n+'Z.193.image_lev1.fits',cache=True) for n in ['02','14','26','38','50']]
    b211_list = [fits.open('/home/mdupont/AIA_data/211_A/aia.lev1_euv_12s.20100613T0539'+n+'Z.211.image_lev1.fits',cache=True) for n in ['02','14','26','38','50']]
    b335_list = [fits.open('/home/mdupont/AIA_data/335_A/aia.lev1_euv_12s.20100613T0539'+n+'Z.335.image_lev1.fits',cache=True) for n in ['02','14','26','38','50']]

    #Block used to store the exposure times for each filter
    expt_171 = []
    expt_193 = []
    expt_211 = []
    expt_335 = []

    #Only care about the CompImageHDU for the fits
    images_171=[fit[1] for fit in b171_list]
    images_193=[fit[1] for fit in b193_list]
    images_211=[fit[1] for fit in b211_list]
    images_335=[fit[1] for fit in b335_list]

    

    #Due to some cards being out of standard format for astropy, we must
    #run a fix function to standardize them before extracting the data
    for image in images_171:
        image.verify('fix')
        expt_171.append(image.header['EXPTIME'])

    for image in images_193:
        image.verify('fix')
        expt_193.append(image.header['EXPTIME'])

    for image in images_211:
        image.verify('fix')
        expt_211.append(image.header['EXPTIME'])

    for image in images_335:
        image.verify('fix')
        expt_335.append(image.header['EXPTIME'])

    #Creates a list of each respective image data from 02s to 50s run time
    concat_171 = [image.data for image in images_171]
    concat_193 = [image.data for image in images_193]
    concat_211 = [image.data for image in images_211]
    concat_335 = [image.data for image in images_335]

    #Will store flux arrays for each filter
    flux_171 = []
    flux_193 = []
    flux_211 = []
    flux_335 = []

    #Will store sliced regions of image to focus on area of interest (The Shock & CME)
    regions_171 = []
    regions_193 = []
    regions_211 = []
    regions_335 = []

    #Must divide the data by its repsective exposure time (Dependent on telescope)
    #in order to get units of DN/s. From there, we define our desired region, cut it out
    #and store it in in list intitilized above
    for i, interest in enumerate(concat_171):
        sf = interest/expt_171[i]
        cutout = Cutout2D(interest,position,size)
        regions_171.append(cutout)
        flux_171.append(sf)

    for i, interest in enumerate(concat_193):
        sf = interest/expt_193[i]
        cutout = Cutout2D(interest,position,size)
        regions_193.append(cutout)
        flux_193.append(sf)

    for i, interest in enumerate(concat_211):
        sf = interest/expt_211[i]
        cutout = Cutout2D(interest,position,size)
        regions_211.append(cutout)
        flux_211.append(sf)

    for i, interest in enumerate(concat_335):
        sf = interest/expt_335[i]
        cutout = Cutout2D(interest,position,size)
        regions_171.append(cutout)
        flux_335.append(sf)
    
    
    #Take the average of the flux for each filter
    #in order to compare the count rates generated by spectrum
    #later 
    avg_flux171 = [np.average(flux) for flux in flux_171]
    avg_flux193 = [np.average(flux) for flux in flux_193]
    avg_flux211 = [np.average(flux) for flux in flux_211]
    avg_flux335 = [np.average(flux) for flux in flux_335]

    


    #Make sure that the region you cutout is indeed your desired region.
    #In this case, we care about the limb where the CME formed
    #for region in regions_193:
        #plt.imshow(region.data, cmap='gray', origin='lower') #stack them
        
    #wave_171[flux] = avg_flux171
    #wave_193[flux] = avg_flux193
    #wave_211[flux] = avg_flux211
    #wave_335[flux] = avg_flux335


    #Maybe graph the average flux values as a function of time
    #from 2s to 50s
    #t= np.linspace(2,50,5)
    #plt.plot(t, avg_flux193)
    #plt.show()


def response_filters():

    #f= open('Emissivity_vs_temp.txt', 'w')

    #Used to store our transmission functions
    trans_171 = []
    trans_193 = []
    trans_211 = []
    trans_335 = []

    
    with open('/home/mdupont/AIA_data/filter_trans/aia_171_ea.txt', 'r') as b_171:
        header = b_171.readline() #Get rid of the header

        

        for line in b_171:
            line = line.strip()
            vals = np.array(line.split(), dtype=float)
            trans_171.append(vals)

    with open('/home/mdupont/AIA_data/filter_trans/aia_193_ea.txt', 'r') as b_193:
        header = b_193.readline()

        

        for line in b_193:
            line = line.strip()
            vals = np.array(line.split(), dtype=float)
            trans_193.append(vals)
    with open('/home/mdupont/AIA_data/filter_trans/aia_211_ea.txt', 'r') as b_211:
        header = b_211.readline()

        

        for line in b_211:
            line = line.strip()
            vals = np.array(line.split(), dtype=float)
            trans_211.append(vals)
    with open('/home/mdupont/AIA_data/filter_trans/aia_335_ea.txt', 'r') as b_335:
        header = b_335.readline()

        

        for line in b_335:
            line = line.strip()
            vals = np.array(line.split(), dtype=float)
            trans_335.append(vals)

    #Turn them into cool numpy arrays instead of simple lists
    trans_171 = np.array(trans_171)
    trans_193 = np.array(trans_193)
    trans_211 = np.array(trans_211)
    trans_335 = np.array(trans_335)

    #Define the effective area and wavelength vectors for each
    #band pass
    ea_171 = (trans_171[:,1])
    wl_171 = trans_171[:,0]

    ea_193 = trans_193[:,1] 
    wl_193 = trans_193[:,0]
    
    ea_211 = trans_211[:,1] 
    wl_211 = trans_211[:,0]

    ea_335 = trans_335[:,1]
    wl_335 = trans_335[:,0]

    LogT = np.linspace(4,9,51)
    Te = 10**LogT

    #only need one set of bins
    wl = {}
    wl['bins'] = wl_171

    #Here are our desired Iron ion stages
    ions = [9,12,24,14,16]
    all_ions = np.arange(27)
    #The indicies corresponde to temperatue hierachy from 2 being 10$^4$
    #and 35 being 10$^7$ in units of Kelvin
    indices = np.arange(2,53)
    all_filters = np.vstack((ea_171,ea_193,ea_211,ea_335))
    aia_responses = ['171','193','211','335']

    filters = {}
    for i, resp in enumerate(aia_responses):
        filt = all_filters[i]
        filters[resp] = filt

    wl_keV = 12.398425/wl_171 #convert to keV
    wl_keV = wl_keV[::-1] #flip it around so it's monotnonically increasing 
    
    ldat = pyatomdb.pyfits.open('/home/mdupont/atomdb/apec_nei_line.fits')
    cdat = pyatomdb.pyfits.open('/home/mdupont/atomdb/apec_nei_comp.fits')
    #This block returns spectra for each iron ion stage as well as temperature and store
    #them in the lists above.
    emp = {}

    spec = {}
    for z1 in all_ions:
        spec[z1] = {}
        emp[z1] = {}
        for k, v in filters.items():
            spec[z1][k] = {}
            for j, index in enumerate(indices):
                ion_spec = pyatomdb.spectrum.make_ion_spectrum(wl_keV,index,Z=26,z1=z1,binunits='keV',
                                                        linefile=ldat,cocofile=cdat,
                                                        dummyfirst=True, nei=True)
                ion_spec = ion_spec[::-1] #Flip the spectrum in order to be correctly folded through the response functions 
                emp[z1][Te[j]] = ion_spec
                mspec = ion_spec*v
                    
                    
                spec[z1][k][Te[j]] = mspec

        
    #Place everything back in Angstroms
    #return spec, wl, filters, emp

    tot_e = {}
    for z1 in all_ions:
        tot_e[z1] = {}
        for k, v in filters.items():
            a = [np.sum(val) for key, val in sorted(spec[z1][k].items())]
            a = np.array(a)
            tot_e[z1][k] = a

    
    #Run test plot of emissvity as function of temperature
    
    ax1 = plt.subplot(221)
    plt.semilogx(Te, tot_e[9]['171'], label='Fe IX')
    plt.semilogx(Te, tot_e[12]['171'], '-+', label='Fe XII')
    plt.semilogx(Te, tot_e[24]['171'],'-x', label='Fe XXIV')
    plt.semilogx(Te, tot_e[14]['171'],'-.', label='Fe XIV')
    plt.semilogx(Te, tot_e[16]['171'],'--', label='Fe XVI')
    plt.title('171 $\AA$ Band')
    

    ax12 = plt.subplot(222)
    plt.semilogx(Te, tot_e[9]['193'], label='Fe IX')
    plt.semilogx(Te, tot_e[12]['193'], '-+', label='Fe XII')
    plt.semilogx(Te, tot_e[24]['193'],'-x', label='Fe XXIV')
    plt.semilogx(Te, tot_e[14]['193'],'-.', label='Fe XIV')
    plt.semilogx(Te, tot_e[16]['193'],'--', label='Fe XVI')
    plt.title('193 $\AA$ Band')

    ax3 = plt.subplot(223)
    plt.semilogx(Te, tot_e[9]['211'], label='Fe IX')
    plt.semilogx(Te, tot_e[12]['211'], '-+', label='Fe XII')
    plt.semilogx(Te, tot_e[24]['211'],'-x', label='Fe XXIV')
    plt.semilogx(Te, tot_e[14]['211'],'-.', label='Fe XIV')
    plt.semilogx(Te, tot_e[16]['211'],'--', label='Fe XVI')
    plt.title('211 $\AA$ Band')

    ax4 = plt.subplot(224)
    plt.semilogx(Te, tot_e[9]['335'], label='Fe IX')
    plt.semilogx(Te, tot_e[12]['335'], '-+', label='Fe XII')
    plt.semilogx(Te, tot_e[24]['335'],'-x', label='Fe XXIV')
    plt.semilogx(Te, tot_e[14]['335'],'-.', label='Fe XIV')
    plt.semilogx(Te, tot_e[16]['335'],'--', label='Fe XVI')
    plt.title('335 $\AA$ Band')
    

    plt.legend()
    plt.show()

    return tot_e
    
    #This block just plots the reponse functions together
    #plt.plot(wl_193, ea_193, label='193')
    #plt.plot(wl_171, ea_171, label='171')
    #plt.plot(wl_211, ea_211, label='211')
    #plt.plot(wl_335, ea_335, label='335x25')
    #plt.title('Filter Response')
    #plt.xlabel('Wavelength in A')
    #plt.ylabel('Effective Area cm$^2$')

    #plt.legend()
    #plt.show()



def multi():
    filters = ['171','193','211','335']
    ions = [9,12,24,14,16]
    Te = np.logspace(4,9,51)
    fe_ions = np.arange(27)

    emis = {}
    with open('total_171.txt','r') as f:
        header=f.readline()
        for line in f:
            i = int(line.split()[1])
            emis[i] = {}
            emis[i]['171'] = np.array(line.split()[2:], dtype=float)
    
    with open('total_193.txt','r') as f:
        header=f.readline()
        for line in f:
            i = int(line.split()[1])
            emis[i]['193'] = np.array(line.split()[2:], dtype=float)
    
    with open('total_211.txt','r') as f:
        header=f.readline()
        for line in f:
            i = int(line.split()[1])
            emis[i]['211'] = np.array(line.split()[2:],dtype=float)
    
    with open('total_335.txt','r') as f:
        header=f.readline()
        for line in f:
            i = int(line.split()[1])
            emis[i]['335'] = np.array(line.split()[2:], dtype=float)


    #grab the emission curves created by response_filters() function

    #Grab ion balances of Fe at different temperatures
    ionbal = []
    for temp in Te:
        bal = pyatomdb.apec.solve_ionbal_eigen(26,temp,False)
        ionbal.append(bal)

    ionbal = np.array(ionbal)

    #Now create dictionary storing the ion_fractions as functions of temperature
    ionf = {}

    for z1 in fe_ions:
        ionf[z1] = {}
        for j, temp in enumerate(Te):
            ionf[z1][temp] = ionbal[j][z1]

    i171 = np.zeros(51, dtype=float)
    i193 = np.zeros(51, dtype=float)
    i211 = np.zeros(51, dtype=float)
    i335 = np.zeros(51, dtype=float)
   

    #171 Filter
    for zed in fe_ions:
      i171 += ionbal[:,zed] * emis[zed]['171']

    for zed in fe_ions:
      i193 += ionbal[:,zed] * emis[zed]['193']

    for zed in fe_ions:
      i211 += ionbal[:,zed] * emis[zed]['211']

    for zed in fe_ions:
      i335 += ionbal[:,zed] * emis[zed]['335']

    


    
    intensity_filters = np.vstack((i171,i193,i211,i335))

    #Create intensity dictionary to pull from later
    intensities = {}
    # intensities['171']=i171

    for k, response in enumerate(filters):
        intensities[response] = intensity_filters[k]


    return emis,ionf,intensities
    
def func(x,Te,c):
    Te = np.logspace(4,9,51)
    emission, ionfrac,it = multi()

    return c*it['171']

def myfit():
    popt,pcov = curve_fit(func,)
      
#Was using this to wrap my head around using the pyatomdb.spectrum.Session(). Just toy spectrum
def toy_spectra():

    ebins = np.linspace(1,2,1001)
    
    a = pyatomdb.spectrum.make_ion_spectrum(ebins,20,Z=26,z1=9, binunits='A')

    Te = np.logspace(4,9, 1000)

    plt.semilogx(Te, a)
    plt.show()




   

    
     