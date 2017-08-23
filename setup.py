from distutils.core import setup, Extension

module1 = Extension('synch',
                    include_dirs = ['/usr/local/include'],
                    libraries = ['m','gsl','gslcblas'],
                    sources = ['synchmodule.c'])

setup (name = 'synchrotron',
       version = '0.01',
       description = 'Synchrotron spectra in python',
       ext_modules = [module1])
