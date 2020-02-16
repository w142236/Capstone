# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:27:11 2020

@author: willk
"""

#IVT = 1/g * Integral_1000^300 (q*V) dP

import numpy as np
from numpy import trapz
from scipy.integrate import simps 
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
#Task 1: Read in the data file 
filename = "C:/Users/willk/Desktop/Capstone.nc"

fh = Dataset(filename, mode =  "r")
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
levs = fh.variables['lev'][:]
time = fh.variables['time'][:] #360 day calendar [1 1 1 0] #47 days (exclusively)
Temp = fh.variables['ta'][:] #47 days (48, 32, 64)
surf_T = fh.variables['tas'][:] #time, lat, lon #Also 2m temp
u_wind = fh.variables['ua'][:] #[48][6][32][64]
v_wind = fh.variables['va'][:]
spec_humidity = fh.variables['hus'][:]#same dims
heights = fh.variables['zg'][:]#same dims
rel_humidity = fh.variables['hur'][:]#same dims

#Obtain sample values for some space and make it a vertical column of air extending from 1000 to 200 hPa

sample_V = np.array([u_wind[12][i][16][32] for i in range(len(levs))], )
sample_q = np.array([spec_humidity[12][i][16][32] for i in range(len(levs))])


## Curves to integrate. Could probably just find the area underneath these curves


# Velocity curve
plt.plot(sample_V,levs)
plt.gca().invert_yaxis()
plt.show()

## Now lets get the area underneath

Area_V = trapz(sample_V, dx = -800) #where dx is from 1000 to 200

print("Area underneath curve is: {}".format(np.round(Area_V)))

#Specific Humidity curve
plt.plot(sample_q,levs)
plt.gca().invert_yaxis()
plt.show()

## Now lets get the area underneath

Area_q = np.abs(trapz(sample_q, dx = -800)) #IMPORTANT: dx is being multiplied twice if we combine these values in the integral(q*V)dP

print("Area underneath curve is: {}".format(np.round(Area_q)))

## Now lets get the sample IVT value

sample_IVT = 1/9.8 * Area_q * Area_V / 800

print("\n\nSo, the IVT value for this vertical column assumed to be 1 x 1 m^2 is (1/g)(Area_q * Area_V)/dx to account for multiplying twice")
print("For this sample the IVT value is: {} kg/ms".format(sample_IVT))