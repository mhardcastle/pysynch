from distutils.core import setup, Extension

module1 = Extension('synch',
                    include_dirs = ['/usr/local/include'],
                    libraries = ['m','gsl','gslcblas'],
                    sources = ['synchmodule3.c'])

setup (name = 'synchrotron',
       version = '0.02',
       description = 'Synchrotron spectra in python3',
       author = 'Martin Hardcastle',
       author_email = 'mjh@extragalactic.info',
       install_requires = [
           'astropy',
           'numpy',
       ],
       ext_modules = [module1],
       py_modules=['synchro'])
