import numpy as np

import math
import matplotlib
import matplotlib.pyplot as plt

from astropy.cosmology import WMAP9
from sncosmo import Source

import astropy.units as u

cosmo=WMAP9
matplotlib.use('WebAgg')

# Model Rabinak&Waxman 2011 doi:10.1088/0004-637X/728/1/63
# ------------------------------------------
# parameters of the model
# * parameters can be modified
# ------------------------------------------

f_p = (0.079+0.13)*0.5
n = 10 # number of solar maxes *
M_sun = 1.989e+30 # [kg]
M = M_sun*n # [kg]
E = 1e+51 # [erg] *
E_51 = E/(1e+51)
R = 1e+13 # [cm] *
R_sun_13 = R/(1e+13)
k = 0.34
k_034 = k/(0.34)
sigma = 5.67e-5 # [erg cm-2 s-1 K-4]
from_Mpc_to_cm = 3.086e+24
from_eV_to_K = 11604.525
h_tagliato = 6.582119569e-16 # [ev s]
h = 2*3.14*h_tagliato*from_eV_to_K # [K s]
c = 2.99792458e+10 # cm/s

# ------------------------------------------

def fluxDensity(t, nu, *args, **kwargs):
    
    argsDict = parseArgs(args, kwargs)
    z = argsDict.pop('z') if 'z' in argsDict else 0.0
    D = cosmo.luminosity_distance(z)*from_Mpc_to_cm # cm
   
    t = (t*86400) # conversion from days to s because the time in chockedclass is given in days
    t = t/(1e+5)
 
    # ------------------------------------------
    # Radius of photosphere
    
    r_ph=3.3e+14*pow(f_p, -0.062)*((pow(E_51,0.41)*pow(k_034,0.093))/(pow((M/M_sun),0.31)))*pow(t,0.81);
    
    # ------------------------------------------
    # Effective Temperature of photosphere
    
    T_ph=(1.6*pow(f_p,-0.037)*((pow(E_51,0.027)*pow(R_sun_13,(1./4.)))/(pow((M/M_sun),0.054)*pow(k_034,0.28)))*pow(t,-0.45))* from_eV_to_K # K
 
    
    # ------------------------------------------
    # Color temperature
    
    T_col = 1.2*T_ph # K
    
    nu = nu*1e-8 # nu converted from A to cm because in chockedclass is given in A
    
    x_lambda = (h*c)/(nu*T_col) # adimensional variable
    
    g_BB = (15./pow(math.pi,4))*(pow(x_lambda,5)/(np.exp(x_lambda)-1.))
    
    f_nu_t = pow((r_ph/D),2)*sigma*pow(T_ph,4)*(T_col/(h*c))*g_BB; # [erg cm-3 s-1]
    
    return f_nu_t



def parseArgs(args, kwargs):
    """
    Parse the arguments to fluxDensity() or intensity(). Supports both
    positional and keyword arguments for now.
    """

    argsDict = kwargs.copy()

    # If there were no extra positional arguments, things are easy.
    if len(args) == 0:
        return argsDict

    
