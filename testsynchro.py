from synchro import SynchSource
from astropy.cosmology import WMAP9
import numpy as np
import matplotlib.pyplot as plt

#s=SynchSource(type='sphere',asph=5,z=0.1,cosmology=WMAP9,spectrum='broken',dpow=1,gbreak=1e3,gmin=10,gmax=1e5,injection=2.1,B=1e-9)
s=SynchSource(type='sphere',asph=5,z=0.1,cosmology=WMAP9,spectrum='aged',ageb=1e-9,age=3e7*365*86400,gmin=50,gmax=1e5,injection=2.0,B=1e-9)

print s.volume

freqs=np.logspace(3,11,40)
emiss=s.emiss(freqs)
plt.plot(freqs,emiss)
plt.xscale('log')
plt.yscale('log')
plt.show()
