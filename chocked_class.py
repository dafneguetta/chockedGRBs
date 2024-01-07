import numpy as np
import math
import chockedday
from chockedday import *
from sncosmo import Source

from scipy.interpolate import (
    RectBivariateSpline as Spline2d
)

from astropy.cosmology import WMAP9
import astropy.units as u

from copy import copy as cp
from math import ceil
from textwrap import dedent

import numpy as np
from astropy import (cosmology, units as u)
from astropy.utils.misc import isiterable
from scipy.interpolate import (
    InterpolatedUnivariateSpline as Spline1d,
    RectBivariateSpline as Spline2d
)

cosmo=WMAP9
def copy(**kwargs):
    assert len(kwargs.values())==1
    return list(kwargs.values())[0]

class ChockedSource(Source):

    _param_names=['input_redshift']
    param_names_latex = ['redshift']
    def __init__(self, phase, wave, redshift=1,
                 time_spline_degree=3, name=None):
        

        self._phase=phase # Assume in days
        self.name = name
        self._wave=wave # Assume wave in Angstroms
        self._parameters=np.array([redshift])
        #self._parameters = np.array([1.])



    def _flux(self, phase, wave):
        """Returns 2-D array of fluxes corresponding to input 1-D arrays of 
        phase and wavelengths.  Phase needs to be in units of days, wave in angstroms"""

        # Get params from property
        input_params={self._param_names[i]: self._parameters[i] for i in range(len(self._param_names))}
        # Get additional params for afterglowpy, as well as convert to recognized names
        #input_params['z']=self._parameters[1]
        input_params['z']=input_params['input_redshift']
        
        # Clean parameters that are not an input to afterglowpy
        #input_params.pop('t0')
        # A check to ensure everyting is correct
        print(input_params)

        # Initialize a time and wavelength array to calculate fluxes over
        t = self._phase # days
        angstrom_array=self._wave
        wave_array=angstrom_array # Angstrom

        fluxes=np.zeros((len(t),len(wave_array)))
        for i,time in enumerate(t):
            # If time is negative, return zero flux
            if time<0:
                fluxes[i]=np.zeros(shape=wave_array.shape)
            else:
                time_array = np.empty(wave_array.shape)
                time_array[:]=time
                
                Fnu = chockedday.fluxDensity(time_array, wave_array, **input_params)
                
                fluxes[i]=Fnu*1e-8 # conversion from erg/s/cm^3 to erg/s/cm^2/A
                
 # Interpolate over fluxes
        f = Spline2d(t, angstrom_array, fluxes, kx=3, ky=3)
       
        self._model_flux = f
        # Retrieve fluxes by inputting the desired phase, wavelength values into interpolation
        #return Fnu(time_array, wave_array)
        return f(phase, wave)

