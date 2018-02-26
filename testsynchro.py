from synchro import SynchSource
from astropy.cosmology import WMAP9
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import matplotlib.pyplot as plt

c=FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
z=0.1
normfreq=1.4e9 # Hz
normflux=1.0 # Hz

s=SynchSource(type='sphere',asph=5,z=z,cosmology=WMAP9,spectrum='broken',dpow=1,gbreak=1e3,gmin=10,gmax=1e5,injection=3.1,B=1e-9,verbose=True)
#s=SynchSource(type='cylinder',arcyl=10,alcyl=10,z=0.1,cosmology=c,spectrum='powerlaw',gmin=10,gmax=1e6,injection=2.0,B=1e-9,verbose=True)

print 'Volume of the object is',s.volume,'m^3'

freqs=np.logspace(3,12,40)
plt.xscale('log')
plt.yscale('log')

print '------- fixed ----------'
s.normalize(normfreq,normflux,method='fixed',bfield=2e-9)
emiss=s.emiss(freqs)
plt.plot(freqs,emiss,label='fixed')
print '------- min energy ----------'
s.normalize(normfreq,normflux,method='minimum_energy',brange=(1e-10,1e-7))
emiss=s.emiss(freqs)
plt.plot(freqs,emiss,label='min_energy')
print '------- equipartition -------'
s.normalize(normfreq,normflux,method='equipartition',brange=(1e-10,1e-7))
emiss=s.emiss(freqs)
plt.plot(freqs,emiss,label='equipartition')


plt.scatter(normfreq*(1+z),normflux*1e-26*s.fnorm/s.volume/s.dfactor,color='red',label='data')
plt.xlabel('Rest-frame frequency (Hz)')
plt.ylabel('Volume emissivity (W m**-3 Hz**-1)')
plt.legend()
plt.show()

