import pyatomdb, numpy

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
  Zlist = range(1,29)
  
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
  for Z in Zlist:
    for z1 in range(1, Z+2):
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
      print(s)
  f.close()
lorentz_power('3.0.8')