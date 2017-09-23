from sherpa.astro.ui import *
import numpy, pyatomdb
from lmfit import Model
import matplotlib.pyplot as plt

# Reading in a bunch of data
intens_tmp = numpy.loadtxt('intensities.txt')
intens={}
intens['171'] = intens_tmp[:,0]
intens['193'] = intens_tmp[:,1]
intens['211'] = intens_tmp[:,2]
intens['335'] = intens_tmp[:,3]
intens['Te'] = numpy.logspace(4,9,51)

# store logs for convenience
intens['log171'] = numpy.log(intens['171'])
intens['log193'] = numpy.log(intens['193'])
intens['log211'] = numpy.log(intens['211'])
intens['log335'] = numpy.log(intens['335'])
intens['logTe']  = numpy.log(intens['Te'] )


# by ion
ion_intens={}
tmp = numpy.loadtxt('total_171.txt', skiprows=1)
ion_intens['171']=tmp[:,2:]

tmp = numpy.loadtxt('total_193.txt', skiprows=1)
ion_intens['193']=tmp[:,2:]

tmp = numpy.loadtxt('total_211.txt', skiprows=1)
ion_intens['211']=tmp[:,2:]

tmp = numpy.loadtxt('total_335.txt', skiprows=1)
ion_intens['335']=tmp[:,2:]
ion_intens['te'] = numpy.log10(numpy.logspace(4,9,51))

# read in the intensities

load_data(1, 'flx171_full.txt',2)
load_data(2, 'flx193_full.txt',2)
load_data(3, 'flx211_full.txt',2)
load_data(4, 'flx335_full.txt',2)


dcache={}

# define my own fit function.
def myfunc(pars, x):
  # This fit function assumes a shock hits at t=t_shock,
  # changing the electron temperature from Te_init to 
  # Te_final
  # density is left as a free parameter.
  # in this case, x is the time
  #  filt is the filter, as an integer (171, 193, 211 or335)
  Te_init = pars[0]
  C = pars[1]
  A = pars[2]
  Te_final = pars[3]
  t_shock = pars[4]
  density = pars[5]
  filt = pars[6]
  print pars
  # get the times before the t_shock
  ret = numpy.zeros(len(x))



  # emissivity of each ion at initial temperature
  emiss_ti = numpy.zeros(27)
  
  for iz in range(27):
    emiss_ti[iz] = 10**numpy.interp(numpy.log(Te_init), \
                                   ion_intens['te'], \
                                   ion_intens['%i'%(filt)][iz])

  # emissivity of each ion at final temperature
  emiss_tf = numpy.zeros(27)
  
  for iz in range(27):
    emiss_tf[iz] = 10**numpy.interp(numpy.log(Te_final), \
                                   ion_intens['te'], \
                                   ion_intens['%i'%(filt)][iz])

  
  # get ionizatoin balance at each temperature, multiply with emissivity
  # to get emission in eacn filter                                 
  itlist = numpy.where(x < t_shock)[0]
  ionbal = pyatomdb.apec.solve_ionbal_eigen(26, Te_init, datacache=dcache)
  for it in itlist:
    ret[it] = sum(emiss_ti * ionbal)
    
  itlist = numpy.where(x >= t_shock)[0]
  
  taulist =  density * (x-t_shock)
  ionbal = pyatomdb.apec.solve_ionbal_eigen(26, Te_final, Te_init=Te_init,\
                                              tau = taulist,\
                                              datacache=dcache)
  
  for it in itlist:
    ret[it] = sum(emiss_tf * ionbal[it,:])

  ret=( A+C*ret)
  return ret

load_user_model(myfunc, 'mfunc1')
add_user_pars("mfunc1", ['Te_init','C','A','Te_final','t_shock','density', 'filt'])

load_user_model(myfunc, 'mfunc2')
add_user_pars("mfunc2", ['Te_init','C','A','Te_final','t_shock','density', 'filt'])


load_user_model(myfunc, 'mfunc3')
add_user_pars("mfunc3", ['Te_init','C','A','Te_final','t_shock','density', 'filt'])

load_user_model(myfunc, 'mfunc4')
add_user_pars("mfunc4", ['Te_init','C','A','Te_final','t_shock','density', 'filt'])

set_source(1,mfunc1)
set_source(2,mfunc2)
set_source(3,mfunc3)
set_source(4,mfunc4)


# set up initial parameters
mfunc1.Te_init=1.8e6
mfunc1.C=4e14
mfunc1.A=25.0
mfunc1.Te_final=4e6
mfunc1.t_shock=4e14
mfunc1.density=1e7

mfunc2.Te_init=mfunc1.Te_init
mfunc2.C=mfunc1.C
mfunc2.A=mfunc1.A
mfunc2.Te_final=mfunc1.Te_final
mfunc2.t_shock=mfunc1.t_shock
mfunc2.density=mfunc1.density
mfunc2.filt=193
mfunc2.filt.freeze()


mfunc3.Te_init=mfunc1.Te_init
mfunc3.C=mfunc1.C
mfunc3.A=mfunc1.A
mfunc3.Te_final=mfunc1.Te_final
mfunc3.t_shock=mfunc1.t_shock
mfunc3.density=mfunc1.density
mfunc3.filt=211
mfunc3.filt.freeze()


mfunc4.Te_init=mfunc1.Te_init
mfunc4.C=mfunc1.C
mfunc4.A=mfunc1.A
mfunc4.Te_final=mfunc1.Te_final
mfunc4.t_shock=mfunc1.t_shock
mfunc4.density=mfunc1.density
mfunc4.filt=335
mfunc4.filt.freeze()



# set some parameter limits
mfunc1.C.min=0.0

mfunc1.Te_init.min=1e4
mfunc1.Te_init.max=1e9


mfunc1.Te_final.min=1e4
mfunc1.Te_final.max=1e9

mfunc1.density.min=1e4

mfunc1.t_shock=900.

mfunc1.t_shock.min=800.
mfunc1.t_shock.max=1000.


mfunc1.filt=171
mfunc1.filt.freeze()

# the following commands make things happen:


fit()

#(does the fit)



#plot_fit(1)
#zzz=raw_input('press a key to display next filter')
#plot_fit(2)
#zzz=raw_input('press a key to display next filter')
#plot_fit(3)
#zzz=raw_input('press a key to display next filter')
#plot_fit(4)
#zzz=raw_input('press a key to exit')
