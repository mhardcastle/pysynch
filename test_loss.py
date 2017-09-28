#!/usr/bin/python

import synch
import numpy as np
import matplotlib.pyplot as plt

values=80
b=0.6e-9
z=0

synch.setspectrum(100,1e6,2.1)
age=np.linspace(0,100,101)
loss=np.zeros_like(age)
for j in range(len(age)):
    synch.setage(age[j]*1e6*365*86400,b)
    loss[j]=synch.loss(1.0,b)

plt.yscale('log')
plt.plot(age,loss)
plt.show()

