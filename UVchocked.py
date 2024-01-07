import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import skysurvey
from shapely import geometry
from skysurvey.tools import utils
from skysurvey.tools import blackbody
from skysurvey.target import Transient
from skysurvey.tools import random_radec
import sncosmo
from astropy.cosmology import WMAP9
import astropy.units as u
from astropy.time import Time
from chocked_class import *
from plotting_lc import my_target_lightcurve
from plotting_lc import my_mag_lightcurve
from plotting_lc import my_show_lightcurve
#size = 5_000 # number of pointings
matplotlib.use('WebAgg')

# Default cosmology used by skysurvey to convert redshift to lum distance
cosmo=WMAP9

# Below are all the requirements for observation with ULTRASAT
# For now, requires ultrasat survey data files to be in same directory as this program
ultrasat_datadir=os.path.dirname(os.path.abspath(__file__))

# Define bandpass for ULTRASAT
wavelength=[]; transmission=[]
# For now, for a given wavelength, get the best transmission value from the 25 radial positions
#with open(ultrasat_datadir+'/transmission.dat' , 'r') as file:
with open('transmission.dat' , 'r') as file:
    transmission_data=file.readlines()
for line in transmission_data:
    values=np.array(line.split(sep=',')).astype(float)
    transmission.append(values[10]) # Get transmission value for tenth radial position

with open(ultrasat_datadir+'/wavelength.dat') as file:
    while True:
        line=file.readline()
        if not line:
            break
        wavelength.append(float(line))
        
# Register band under name 'ultrasat'
band=sncosmo.Bandpass(wavelength, transmission, name='ultrasat')
sncosmo.register(band)

# Get zp and skynoise for ultrasat
with open(ultrasat_datadir+'/zp.dat') as file:
    ultrasat_zp=np.array(file.readline().split(sep=',')).astype(float)[10]
with open(ultrasat_datadir+'/varperpix.dat') as file:
    ultrasat_noise=np.array(file.readline().split(sep=',')).astype(float)[10] # This is variance of electron count per pixel

# Both currently use the tenth radial position given

# Need negative values for t so we properly generate zero fluxes at negative phase.
t = np.concatenate((np.array([-1.0e4, -10, -1]), np.geomspace(1.0e3, 1.0e7, 100)),axis=0) # in seconds

#t = np.geomspace(1.0e3, 1.0e7, 100)
# Needs to encompass entire range ultrasat band is defined on
#angstrom_array=np.concatenate((np.linspace(1000, 2290, num=100),
 #                              np.linspace(2300, 2900, num=250),
  #                             np.linspace(2910, 11100, num=100)), axis=0) # In Angstroms
angstrom_array= np.linspace(1000, 11100,100000)
phases=t/86400 # Convert to days

model=dict(redshift={"kwargs":{"zmax":0.5}, "as":"z"},
           input_redshift={"func": copy,
                "kwargs":{"a":"@z"}},
            t0 = {"func": np.random.uniform,
                 "kwargs": {"low":60_000, "high":61_200}},
            radec={"func":random_radec,
                  "as":["ra","dec"]})



#Define a base Transient class
class ChockTransient(skysurvey.TSTransient):
    _KIND = "chock_test" # your transient name
    _TEMPLATE = sncosmo.Model(ChockedSource(phases, angstrom_array))# skysurvey.Template or sncosmo.source (or name)
    _RATE = 1000 # float (Volume Rate in Gyr/year ; or redshift function)
    _MODEL = model

my_transient=ChockTransient()

# Define ULTRASAT parameters
footprint = geometry.Point(0,0).buffer(8.05823906207) # Assume circular FOV with area of 204 deg^2

time_start=Time('2026-06-01')
time_start.format='mjd'
mjd_range = [time_start.value, time_start.value+30] # First date is June 01, 2026, 1 month range

data = {}
steps=np.round((mjd_range[1]-mjd_range[0])/0.00347222222)
data["mjd"] = np.linspace(mjd_range[0], mjd_range[0]+0.00347222222*steps, num=int(np.round((mjd_range[1]-mjd_range[0])/0.00347222222))) # Latter is 5 min interval in MJD units
npointings = len(data["mjd"])
data["ra"] = np.array([220]*npointings)
data["dec"] = np.array([66]*npointings) # Uses N1 field, change to 42, -66 for S1
data["skynoise"] = np.random.normal(loc=ultrasat_noise, scale=0, size=npointings) # Currently zero scatter
data["gain"] = 2
data["zp"] = ultrasat_zp
data["band"] = "ultrasat"


mysurvey = skysurvey.Survey.from_pointings(data, footprint=footprint)


# Set skyarea for ultrasat survey
ultrasat_area=geometry.Point(220,66).buffer(8.06489662848)
#mysurvey.show()
# Generate target data

result=my_transient.from_draw(size=20,
                             tstart="2026-06-01", tstop="2026-07-01",
                             skyarea=ultrasat_area) # One month range



print(result.data)
print(mysurvey.data)
dset = skysurvey.DataSet.from_targets_and_survey(result, mysurvey)
#fig = dset.show_target_lightcurve(phase_window=[-1,10], show_truth=True)
#fig2 = dset.show_lightcurve(band="ultrasat", in_mag=True)
my_target_lightcurve(dset, phase_window=[-1,10])
my_show_lightcurve(result, band="ultrasat", index=3, params=None,
                            ax=None, fig=None, colors=None,
                            time_range=[-1,10], npoints=500,
                            zp=None, zpsys="ab",
                            format_time=True, t0_format="mjd", 
                            in_mag=False)
plt.ylim(1.0e-6,1.0e-1)
#my_mag_lightcurve(result, band="ultrasat", index=10, in_mag=True)
my_mag_lightcurve(result, band="ultrasat", index=3, params=None,
                            ax=None, fig=None, colors=None,
                            time_range=[-1,10], npoints=500,
                            format_time=True, t0_format="mjd", 
                            in_mag=True)


plt.ylim((25,17))
#print(dset.data)
plt.show()
