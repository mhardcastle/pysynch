from __future__ import division
from __future__ import print_function
from builtins import object
import numpy as np

PI=np.pi
PARSEC=3.08568e16
JANSKY=1.0e-26
V_C=299792458.0
MU_0=4.0e-7*PI

import synch

class SynchSource(object):
    '''Class that sets up a physical object and normalizes a synchrotron
    spectrum based on its observable properties. Methods generally
    update attributes, see method docstrings for details.
    '''

    def _init_distances(self):
        # scale in kpc/arcsec
        if self.scale is None:
            self.scale=1.0/(self.cosm.arcsec_per_kpc_proper(self.z).value)
        # flux normalization
        if self.fnorm is None:
            self.fnorm=4*PI*(1e6*PARSEC*self.cosm.luminosity_distance(self.z).value)**2.0/(1.0+self.z)
    
    def arcsec_to_metres(self,theta):
        self._init_distances()
        return theta*self.scale*1000.0*PARSEC
    
    def __init__(self,type='sphere',cosmology=None,z=0,cmbtemp=2.725,verbose=False,gmin=None,gmax=None,injection=None,spectrum='powerlaw',**kwargs):
        if gmin is None or gmax is None or injection is None:
            raise RuntimeError('gmin, gmax, injection must be specified')
        self.verbose=verbose
        # sort out cosmology. This should be an astropy cosmology object
        self.cosm=cosmology
        # these are set up so we can easily apply a fixed
        # non-cosmological distance later, as in synch
        self.scale=None
        self.fnorm=None
        self.z=z
        self.type='type'
        # set up the volume
        if type=='sphere':
            if 'rsph' in kwargs:
                self.rsph=kwargs['rsph']
            else:
                self.rsph=self.arcsec_to_metres(kwargs['asph'])
            self.volume=4*PI*self.rsph**3.0/3.0
        elif type=='cylinder':
            if 'lcyl' in kwargs:
                self.lcyl=kwargs['lcyl']
                sefl.rcyl=kwargs['rcyl']
            else:
                self.lcyl=self.arcsec_to_metres(kwargs['alcyl'])
                self.rcyl=self.arcsec_to_metres(kwargs['arcyl'])
            self.volume=PI*self.rcyl**2.0*self.lcyl
        elif type=='ellipsoid':
            if 'major' in kwargs:
                self.major=kwargs['major']
                self.minor=kwargs['minor']
            else:
                self.major=self.arcsec_to_metres(kwargs['amajor'])
                self.minor=self.arcsec_to_metres(kwargs['aminor'])
            self.volume=PI*self.minor**2.0*self.major/6.0

        else:
            raise NotImplementedError('geometry '+type)

        # beaming -- not yet implemented
        self.doppler=1.0
        self.dfactor=1.0

        # synchrotron spectrum
        self.synchnorm=1.0
        self.gmin=gmin
        self.gmax=gmax
        self.injection=injection
        if 'B' in kwargs:
            self.B=kwargs['B']
        else:
            self.B=None

        self.spectrum=spectrum
        if spectrum=='powerlaw':
            pass
        elif spectrum=='broken':
            self.gbreak=kwargs['gbreak']
            self.dpow=kwargs['dpow']
        elif spectrum=='aged':
            self.age=kwargs['age']
            self.ageb=kwargs['ageb']
        else:
            raise NotImplementedError('spectrum '+spectrum)
            
    def _setsynch(self):
        synch.setspectrum(self.gmin,self.gmax,self.injection)
        if self.spectrum=='powerlaw':
            synch.setpow()
        elif self.spectrum=='broken':
            synch.setbreak(self.gbreak,self.dpow)
        elif self.spectrum=='aged':
            synch.setage(self.age,self.ageb)
        else:
            raise NotImplementedError('spectrum '+spectrum)
            
    def emiss(self,freq):
        self._setsynch()
        freq=np.array(freq)
        it=np.nditer([freq,None])
        for f,r in it:
            r[...]=synch.emiss(self.synchnorm,self.B,f)
        return it.operands[1]
    
    def cmb_ic_emiss(self,freq):
        self._setsynch()
        freq=np.array(freq)
        it=np.nditer([freq,None])
        for f,r in it:
            r[...]=synch.cmb_ic_emiss(self.synchnorm,f,self.z)
        return it.operands[1]

    def normalize(self,freq,flux,zeta=1.0,method=None,tol=1e-6,**kwargs):
        if method is None:
            raise RuntimeError('Method must be specified with method keyword.')
        self._init_distances()
        nu=freq*(1+self.z)/self.doppler
        wem=flux*JANSKY*self.fnorm/self.volume/self.dfactor
        # compute single-electron norm and energy density
        self._setsynch()
        tno=synch.intne(1)
        if self.verbose: print('Total no of electrons (pre-norm) is %g m^-3' % tno)
        eln=1.0/tno
        if self.verbose: print('Total no of electrons (post-norm) is %f m^-3' % synch.intne(eln))
        ed=synch.intene(eln)
        if self.verbose: print('Energy density (one electron) is %g J m^-3' % ed)
        if self.verbose: print('Wanted emission rate at %g Hz: %g W/Hz/m^3' % (nu,wem))
        if method=='fixed':
            bfield=kwargs['bfield']
            bed=bfield**2.0/(2.0*MU_0)
            eo=synch.emiss(eln,bfield,nu)
            norm=(wem/eo)*eln
            self.B=bfield
            self.synchnorm=norm
            self.electron_energy_density=norm*ed/eln
            self.bfield_energy_density=bed
            self.total_energy_density=self.electron_energy_density+self.bfield_energy_density
        elif method=='equipartition':
            bmin,bmax=kwargs['brange']
            if self.verbose: print('Finding the B-field such that %f * B^2/mu_0 = total electron energy' % zeta)

            # check the bounds
            bfield=bmax
            bed=bfield**2.0/(2.0*MU_0)
            norm=eln*zeta*bed/ed

            eu=synch.emiss(norm,bfield,nu)
      
            bfield=bmin
            bed=bfield**2.0/(2.0*MU_0)
            norm=eln*zeta*bed/ed
            el=synch.emiss(norm,bfield,nu)
      
            if self.verbose: print("Emission rate limits are %g -- %g W Hz^-1 m^-3" % (el,eu))

            if wem<el or wem>eu:
                raise RuntimeError('Not bracketing a root; can\'t search')

            emid=0
            while (abs(emid-wem)/wem)>tol:
                bfield=np.exp((np.log(bmin)+np.log(bmax))/2.0)
                bed=bfield**2.0/(2.0*MU_0)
                norm=eln*zeta*bed/ed
                emid=synch.emiss(norm,bfield,nu)
                if wem<emid:
                    bmax=bfield
                else:
                    bmin=bfield

            self.B=bfield
            self.synchnorm=norm
            self.electron_energy_density=norm*ed/eln
            self.bfield_energy_density=bed
            self.total_energy_density=self.electron_energy_density+self.bfield_energy_density
        elif method=='minimum_energy':
            bmin,bmax=kwargs['brange']
            if self.verbose: print('Finding the B-field that minimizes B^2/mu_0  + %f * total electron energy' % zeta)
            bfield=bmin
            bed=bfield**2.0/(2.0*MU_0)
            el=synch.emiss(eln,bfield,nu)
            norm=wem/el
            el=bed+zeta*norm*ed
            if self.verbose: print("Lower B value total energy density is %g J m^-3" % el)

            bfield=bmax
            bed=bfield**2.0/(2.0*MU_0)
            eu=synch.emiss(eln,bfield,nu)
            norm=wem/eu
            eu=bed+zeta*norm*ed
            if self.verbose: print("Upper B value total energy density is %g J m^-3" % eu)

            # pick an initial central point using the golden ratio -- in log space 

            GR=0.61803399
            GC=(1.0-GR)
            TOL=1.0e-6
            bfield=np.exp(np.log(bmin)+GC*(np.log(bmax/bmin)))
            if self.verbose: print("Initial midpoint is %g T" % bfield)
            bed=bfield**2.0/(2.0*MU_0)
            emid=synch.emiss(eln,bfield,nu)
            norm=wem/emid
            emid=bed+zeta*norm*ed
            if self.verbose: print("Mid-point total energy density is %g J m^-3" % emid)
            if emid>eu or emid>el:
                raise RuntimeError("Not bracketing a minimum; can't search.")

            b0=np.log(bmin)
            b3=np.log(bmax)

            b1=np.log(bfield)
            b2=b1+GC*(b3-b1)
            f1=emid
            bfield=np.exp(b2)
            bed=bfield**2.0/(2.0*MU_0)
            emid=synch.emiss(eln,bfield,nu)
            norm=wem/emid
            f2=bed+zeta*norm*ed
            while (abs(b3-b0)>TOL*(abs(b1)+abs(b2))):
                if (f2<f1):
                    b0=b1
                    b1=b2
                    b2=GR*b1+GC*b3
                    f1=f2
                    bfield=np.exp(b2)
                    bed=bfield**2.0/(2.0*MU_0)
                    emid=synch.emiss(eln,bfield,nu)
                    norm=wem/emid
                    f2=bed+zeta*norm*ed
                else:
                    b3=b2
                    b2=b1
                    b1=GR*b2+GC*b0
                    f2=f1
                    bfield=np.exp(b1)
                    bed=bfield**2.0/(2.0*MU_0)
                    emid=synch.emiss(eln,bfield,nu)
                    norm=wem/emid
                    f1=bed+zeta*norm*ed

            if (f1<f2):
                bfield=np.exp(b1)
            else:
                bfield=np.exp(b2)
            if self.verbose: print("Minimum-energy B-field is %g T" % bfield)
            bed=bfield**2.0/(2.0*MU_0)
            self.B=bfield
            emid=synch.emiss(eln,bfield,nu)
            norm=wem/emid
            self.synchnorm=norm*eln
            self.electron_energy_density=norm*ed
            self.bfield_energy_density=bed
            self.total_energy_density=zeta*self.electron_energy_density+self.bfield_energy_density
        else:
            raise NotImplementedError('method '+method)
            
        if self.verbose:
            print('Field fit %g T' % self.B)
            print('B-field energy density is %g J/m^3' % self.bfield_energy_density)
            print('Normalized total number density of electrons is %g m^-3' % synch.intne(self.synchnorm))
            print('Electron energy density is %g J/m^3' % self.electron_energy_density)
            print('Total energy density is %g J/m^3' % self.total_energy_density)
                
