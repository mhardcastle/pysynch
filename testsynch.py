#!/usr/bin/python

import synch
import numpy as np
import matplotlib.pyplot as plt

values=80
b=0.6e-9

synch.setspectrum(100,1e6,2.1)

freq=np.logspace(3,11,values)
emiss=np.zeros(values)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Emissivity')

for j in range(0,11):

    age=4*j
    synch.setage(age*1e6*365*86400,b)

    for i in range(0,values):
        emiss[i]=synch.emiss(1.0,b,freq[i])
        # print '%g %g' % (freq[i],emiss[i])

    plt.plot(freq,emiss,label=('Age %g Myr' % age))

plt.legend(loc=0)
plt.show()
