from synchro import SynchSource
#from astropy.cosmology import WMAP9 as c
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import matplotlib.pyplot as plt

c=FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

#s=SynchSource(type='sphere',asph=5,z=0.1,cosmology=WMAP9,spectrum='broken',dpow=1,gbreak=1e3,gmin=10,gmax=1e5,injection=2.1,B=1e-9)
s=SynchSource(type='sphere',asph=10,z=0.1,cosmology=c,spectrum='powerlaw',gmin=10,gmax=1e6,injection=2.0,B=1e-9,verbose=True)

print s.volume
print s.fnorm

s.normalize(1.4,1.0,method='equipartition',brange=(1e-10,1e-8))

'''

freqs=np.logspace(3,11,40)
emiss=s.emiss(freqs)
plt.plot(freqs,emiss)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
