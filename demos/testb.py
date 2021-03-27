from __future__ import division
from __future__ import print_function
from past.utils import old_div
from synchro import SynchSource
from astropy.cosmology import WMAP9
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import matplotlib.pyplot as plt

c=FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
z=0.1
normfreq=1.4e9 # Hz
normflux=1.0 # Hz

s=SynchSource(type='sphere',asph=5,z=z,cosmology=WMAP9,spectrum='broken',dpow=1,gbreak=1e3,gmin=10,gmax=1e5,injection=2.5,B=1e-9,verbose=True)
#s=SynchSource(type='cylinder',arcyl=10,alcyl=10,z=0.1,cosmology=c,spectrum='powerlaw',gmin=10,gmax=1e6,injection=2.0,B=1e-9,verbose=True)

print('Volume of the object is',s.volume,'m^3')

s.normalize(normfreq,normflux,method='equipartition',brange=(1e-10,1e-7))

bfield=s.B

bvals=s.B*np.logspace(-1,1,41)
energy=np.zeros_like(bvals)

for i,b in enumerate(bvals):
    
    s.normalize(normfreq,normflux,method='fixed',bfield=b)
    energy[i]=s.total_energy_density

plt.plot(bvals,old_div(energy,np.min(energy)))
plt.xscale('log')
plt.yscale('log')

plt.show()

