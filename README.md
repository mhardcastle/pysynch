# pysynch
Python interface to synchrotron libraries

This provides an interface to the C synchrotron libraries used by
Hardcastle et al (1998) and in subsequent work, and python wrappers
for some basic functions (synchro, see below). pysynch depends on the
integration routines in the GNU Scientific Library (GSL) and synchro
depends on numpy and (indirectly) on astropy as well. Building these
modules has been mostly tested on Linux but it reportedly also works
on Macs if GSL is installed (e.g. through brew).

## installation

To install pysynch clone the git repository and do

```
python setup.py install
```

(the usual options to install may be used).

Then you should be able to

```
import synch
```

Some test code (e.g. `testsynch.py`, `test_ic.py`, `test_loss.py`) is provided in the `demos` directory.

## warning

pysynch is *not* thread-safe.

## pysynch functions

The synch library itself provides the following low-level functions as
interfaces to the C code.

Units are T for magnetic field, Hz for frequency, s for time, and
Lorentz factor for electron energy.

* `setspectrum(gamma_min,gamma_max,power_law_index)`: sets the basic
  electron spectrum parameters. Should be called before anything else:
  there are defaults but they are probably not what you want.

* `setpow()`: sets a basic power-law model. The default; only needed
  if you are changing models.

* `setbreak(gamma_break, delta)`: use a broken power-law electron
  spectrum with a break of delta (that is the amount by which the
  spectrum steepens) at gamma_break.

* `setage(age, agefield)`: use a J-P aged spectrum with age age and
  ageing field agefield. Note that it's your responsibility to make
  sure that ageb takes account of the CMB if that's what you want.

* `emiss(norm, field, frequency)`: calculate the synchrotron
  emissivity J(nu). norm is the normalization of the electron energy
  spectrum, i.e. N_0 in N(E) = N_0 E**-p

* `cmb_ic_emiss(norm, frequency, z)`: calculate the inverse-Compton
  emissivity J_ic(nu). The redshift z must be specified.

* `intne(norm)` and `intene(norm)`: integrate the electron energy
  spectrum over the energy range.

## synchro

This code provides some of the functionality of the synch code of
Hardcastle et al (1998). It allows the user to set up one or more
Python instances which represent real synchrotron/IC sources and carry
out certain operations on them. SSC parts of the code are not yet implemented. 
The code depends on numpy and (indirectly) on astropy.

To use it install the package as above and then

```
from synchro import SynchSource
```

Examples of its use are in `demos/testsynchro.py`.

### creating an instance

To set up a source to model create an instance of the SynchSource
object, e.g.

```
ss=SynchSource(type='sphere', gmin=10, gmax=100000, z=0.1, injection=2.0,
spectrum='powerlaw', cosmology=my_cosmology, asph=10)
```

`cosmology` here is an astropy cosmology instance, or something equivalent, which provides the
methods `arcsec_per_kpc_proper`, `arcsec_to_metres` and
`cosm.luminosity_distance`. The redshift `z` will be used with this instance.

Source types may be 'sphere', 'cylinder' or 'ellipsoid'. Each has its
own parameters for sizes. Either physical or angular sizes may be
used, but *not* a mixture.

* sphere: use keywords `rsph` for the radius of the object in metres
  or `asph` for arcsec.

* cylinder: use `lcyl` and `rcyl` (metres) or `alcyl` and `arcycl`
 (arcsec) for the length and radius.

* ellipsoid: use `major` and `minor` (metres) or `amajor` and `aminor`
  (arcsec).

Spectra may be 'broken', 'powerlaw' or 'aged'. Each again has its own
parameters to pass to SynchSource:

* powerlaw: None

* broken: `gbreak` and `dpow`, the Lorentz factor of the break and the
  change in power-law index

* aged: `age` and `ageb`, the age (in seconds) and ageing field (in
  T), and optionally `cmbage` which must evaluate to `True` if you
  wish to include inverse-Compton losses at the redshift of the source
  in the ageing calculation.

Set the magnetic field in the call to SynchSource with the `B` keyword
(or set the B attribute of the instance).

Set `verbose=True` to get printed reports of what the various methods
do as they happen.

Initialization sets the instance attribute `volume` (in m^3), `scale` (in kpc/arcsec) and `fnorm` (the conversion factor between flux density in W/Hz/m**2 and luminosity in W/Hz).

### instance methods

* `emiss(freq)`: calculate the synchrotron emissivity at a frequency or a list/array of
 frequencies (in the rest frame of the object) and return results. The
 `B` attribute of the instance must be set. `emiss` uses the
 `synchnorm` attribute which will be set to an arbitrary value unless
 `normalize` has been called.

* `cmb_ic_emiss(freq)`: calculate the CMB IC emissivity at a frequency or a list/array of
 frequencies (in the rest frame of the object) and return results. The
 `z` attribute of the instance must be set. `cmb_ic_emiss` uses the
 `synchnorm` attribute which will be set to an arbitrary value unless
 `normalize` has been called.

* `normalize(frequency, flux, method='METHOD')`: set the synchrotron normalization by one of several
  possible methods using an observation of flux density `flux` at
  observer-frame frequency `frequency`.

  Possible normalization methods include:

  * 'fixed': a parameter `bfield` must be passed -- the electron normalization
    is then adjusted to produce the observed emission at that field strength.

  * 'equipartition': Find the equipartition field within a given range.
    A tuple `brange` of the minimum and maximum field
    strengths to use must be passed to the method.

  * 'minimum_energy': Find the minimum-energy field within a given range.
    A tuple `brange` of the minimum and maximum field
    strengths to use must be passed to the method.

  For methods 'equipartition' and 'minimum_energy' the parameter `zeta` is one plus the ratio of energy density in radiating to non-radiating particles (often denoted kappa). The default value is 1.0, i.e. there are no non-radiating particles. Values < 1.0 are not physically meaningful. For true minimum energy `zeta` should always be left at the default value of 1.0.

  normalize sets the instance attributes `B`, `synchnorm`, `bfield_energy_density`,
  `electron_energy_density` and `total_energy_density`, which are
  stored in SI units. `total_energy_density` is by definition equal to `bfield_energy_density + zeta * electron_energy_density`.

## synchro tutorial

Suppose we want to estimate the equipartition magnetic field in the E
lobe of Cygnus A. This source is at a redshift of 0.0565 and we measure a flux of 144.42 Jy at 4.525 GHz. We approximate the shape of the lobe by an ellipse with major axis 50 arcsec and minor axis 30 arcsec. We take the electron spectrum to be a power law with gamma_min=1 and gamma_max=1e5.

```
from synchro import SynchSource
from astropy.cosmology import FlatLambdaCDM
c=FlatLambdaCDM(H0=70, Om0=0.3)

s=SynchSource(type='ellipsoid',gmin=1, gmax=1e5, z=0.0565, injection=2.0, spectrum='powerlaw', cosmology=c, amajor=50, aminor=30)

s.normalize(4.525e9,144.42,method='equipartition',brange=(1e-10,1e-7))

print(s.B)
```

This gives a magnetic field strength estimate of 4.3 nT. 

Now we decide that we prefer to model with a broken power-law electron spectrum where the energy index steepens from 2 to 3:

```
s=SynchSource(type='ellipsoid',gmin=1, gmax=1e5, gbreak=6000, dpow=1, z=0.0565, injection=2.0, spectrum='broken', cosmology=c, amajor=50, aminor=30)

s.normalize(4.525e9,144.42,method='equipartition',brange=(1e-10,1e-7))

print(s.B)
```

We see that this changes (increases) the magnetic field to 4.5 nT -- which is expected since the normalizing frequency is above the break we have put in.

We find the total energy density in the lobes in a similar way

```
print(s.total_energy_density)
```

to get a value of 1.6e-11 J/m**3. The total pressure in Pa is given by dividing by 3:

```
print(s.total_energy_density/3)
```

so the equipartition pressure in the lobes is 5.4e-12 Pa.

We can plot the radio spectral luminosity as a function of (source frame) frequency:

```
import numpy as np
import matplotlib.pyplot as plt

frequencies=np.logspace(6,11,100)
plt.plot(frequencies,s.emiss(frequencies)*s.volume)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Luminosity (W/Hz)')
plt.xlabel('Frequency (Hz)')
plt.show()
```

Now we decide we want to plot the predicted flux density as a function of observer-frame frequency, labelling the normalizing point:

```
plt.plot(frequencies,1e26*s.emiss(frequencies*(1+0.0565))*s.volume/s.fnorm)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Flux (Jy)')
plt.xlabel('Frequency (Hz)')
plt.scatter(4.525e9,144.42,color='red')
plt.show()

```
