#!/usr/bin/python

from builtins import range
import synch
import numpy as np
import matplotlib.pyplot as plt

values=80
b=0.6e-9
z=0

synch.setspectrum(100,1e6,2.1)

synch_freq=np.logspace(3,11,values)
synch_emiss=np.zeros(values)
ic_emiss=np.zeros(values)
ic_freq=2.4e17*np.logspace(-4,4,values)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Emissivity')

for j in range(0,11):

    age=5*j
    synch.setage(age*1e6*365*86400,b)

    for i in range(0,values):
        #synch_emiss[i]=synch.emiss(1.0,b,synch_freq[i])
        ic_emiss[i]=synch.cmb_ic_emiss(1.0,ic_freq[i],z)
        # print '%g %g' % (freq[i],emiss[i])

    #plt.plot(synch_freq,synch_emiss,label=('Age %g Myr' % age))
    plt.plot(ic_freq,ic_emiss,ls='--',label='IC age %g Myr' % age)

plt.xlim(3e17,3e21)
plt.legend(loc=0)
plt.show()
