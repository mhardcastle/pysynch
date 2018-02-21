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
    
    def __init__(self,type='sphere',cosmology=None,z=0,cmbtemp=2.73,verbose=False,gmin=None,gmax=None,injection=None,spectrum='powerlaw',**kwargs):
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
    
    def normalize(self,freq,flux,zeta=1.0,method='equipartition',**kwargs):
        self._init_distances()
        nu=freq*(1+self.z)/self.doppler
        wem=flux*JANSKY*self.fnorm/self.volume/self.dfactor
        # compute single-electron norm and energy density
        self._setsynch()
        tno=synch.intne(1)
        print 'Total no of electrons (pre-norm) is %g m^-3' % tno
        eln=1.0/tno
        print 'Total no of electrons (post-norm) is %f m^-3' % synch.intne(eln)
        ed=synch.intene(eln)
        print 'Energy density (one electron) is %g J m^-3' % ed
        if self.verbose:
            print 'Wanted emission rate at %g Hz: %g W/Hz/m^3' % (nu,wem)
        if method=='equipartition':
            bmin,bmax=kwargs['brange']
            print 'Emissivity (one electron) at bmin is %g' % synch.emiss(eln,bmin,nu)
            if self.verbose: print 'Finding the B-field such that %f * B^2/mu_0 = total electron energy' % zeta

            # check the bounds
            BFIELD=bmax
            bed=BFIELD**2.0/(2.0*MU_0)
            norm=eln*(zeta*bed/ed)
            print bed,ed,eln,norm
            print bed,synch.intene(norm)
            eu=synch.emiss(norm,BFIELD,nu);
      
            BFIELD=bmin;
            bed=BFIELD**2.0/(2.0*MU_0)
            norm=eln*(zeta*bed/ed)
            el=synch.emiss(norm,BFIELD,nu)
            print bed,ed,eln,norm
      
            print "Emission rate limits are %g -- %g W Hz^-1 m^-3" % (el,eu)

            if wem<el or wem>eu:
                raise RuntimeError('Not bracketing a root; can\'t search')
            
            for i in range(20):
                BFIELD=np.exp((np.log(bmin)+np.log(bmax))/2.0)
                bed=BFIELD**2.0/(2.0*MU_0)
                norm=eln*(zeta*bed/ed)
                emid=synch.emiss(norm,BFIELD,nu)
                print bed,ed,eln,norm
                if wem<emid:
                    bmax=BFIELD
                else:
                    bmin=BFIELD

            print 'Field fit %g T' % BFIELD
            
            print 'B-field energy density is %g J/m^3' % bed
            print 'Normalized total no electrons is %g m^-3' % synch.intne(norm)
            print 'Electron energy density is %g J/m^3' % synch.intene(norm)
