import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
from lmfit import Model
from numpy import sqrt, pi, exp



resp = ['171','193','211','335']
# read in the data
intens_tmp = np.loadtxt('/home/mdupont/2017_Projects/ionica/NEI/shocks/intensities.txt')
intens={}
intens['171'] = intens_tmp[:,0]
intens['193'] = intens_tmp[:,1]
intens['211'] = intens_tmp[:,2]
intens['335'] = intens_tmp[:,3]
intens['Te'] = np.logspace(4,9,51)



# store logs for convenience
intens['log171'] = np.log(intens['171'])
intens['log193'] = np.log(intens['193'])
intens['log211'] = np.log(intens['211'])
intens['log335'] = np.log(intens['335'])
intens['logTe']  = np.log(intens['Te'])
signal_tmp = np.loadtxt('/home/mdupont/2017_Projects/ionica/NEI/shocks/allflx.txt')

signal={}
signal['171'] = signal_tmp[:,0]
signal['193'] = signal_tmp[:,1]
signal['211'] = signal_tmp[:,2]
signal['335'] = signal_tmp[:,3]

# made up 20 minute time intervals. Not relevant for now
signal['time'] = np.linspace(0,2*60.0*(len(signal['171'])-1), len(signal['171']))

# dummy array of xs. We don't have a function which varies as x, so can ignore
# this

timmy = np.linspace(0,1788,147,dtype=int)
print signal['time'].size
xv = [1,2,3,4]

def myfitfunc(x, Te, C, a):
  # fit function. INterpolates on the temperature grid to find the predicted
  # emissivity at Te, multiplies by C, and returns the flux in the 4 filters
  # in order [171,193,211,335]
  
  ret = []
  
  for filt in ['171','193','211','335']:
    ret.append(a+C*np.exp(np.interp(np.log(Te), \
                             intens['logTe'], \
                             intens['log%s'%(filt)],\
                             left = 0.0,\
                             right = 0.0)))
  return np.array(ret)

def gaud(x, Te,c,a):
    return Te*x**2+c*x+a
    

gmod = Model(myfitfunc)

for i in range(147):
#print [signal['171'][i], \
#                           signal['193'][i], \
#                           signal['211'][i], \
#                           signal['335'][i]]
  k=sopt.curve_fit(myfitfunc, xv, [signal['171'][i], \
                             signal['193'][i], \
                             signal['211'][i], \
                             signal['335'][i]], p0=[1.8e6,4e14, 5.0e3])

  j = gmod.fit([signal['171'][i], signal['193'][i],signal['211'][i],signal['335'][i]],x=xv,Te=1.8e6,C=4e14,a=5.0e3)

  
  print i, k[0] 

  print j.best_fit

  res['Te'][i] = k[0][0]
  res['norm'][i] = k[0][1]
  res['offset'][i] = k[0][2]
result1 = gmod.fit(signal['171'],x=timmy,Te=1.8e6,c=4e14,a=1e3)
result2 = gmod.fit(signal['193'],x=timmy,Te=1.8e6,c=4e14,a=1e3)
result3 = gmod.fit(signal['211'],x=timmy,Te=1.8e6,c=4e14,a=1e3)
result4 = gmod.fit(signal['335'],x=timmy,Te=1.8e6,c=4e14,a=1e3)


plt.plot(timmy,signal['171'],'c--')
plt.plot(timmy,signal['193'],'g--')
plt.plot(timmy,signal['211'],'r--')
plt.plot(timmy,signal['335'],'b--')

plt.plot(timmy,result1.best_fit,'c')
plt.plot(timmy,result2.best_fit,'g')
plt.plot(timmy,result3.best_fit,'r')
plt.plot(timmy,result4.best_fit,'b')

plt.show()
