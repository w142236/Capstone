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
filename = "C:/Users/willk/Code_Research/Capstone/Capstone.nc"

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

sample_V = np.array([u_wind[0][i][16][33] for i in range(len(levs))], )
sample_q = np.array([spec_humidity[0][i][16][33] for i in range(len(levs))])


## Curves to integrate. Could probably just find the area underneath these curves


# Velocity curve
plt.plot(sample_V,levs)
plt.gca().invert_yaxis()
plt.show()

## Now lets get the area underneath
dP = -800 

Area_V = trapz(sample_V, dx = dP) #where dx is from 1000 to 200
#Area_V_Simpson = simps(sample_)
print("Trapezoidal Area underneath curve is: {}".format(np.round(Area_V)))
#print("Using Simpsons rule, Area underneath the curve is: {}".format(np.round()))
#Specific Humidity curve
plt.plot(sample_q,levs)
plt.gca().invert_yaxis()
plt.show()

## Now lets get the area underneath

Area_q = np.abs(trapz(sample_q, dx = dP)) #IMPORTANT: dx is being multiplied twice if we combine these values in the integral(q*V)dP

print("Trapezoidal Area underneath curve is: {}".format(np.round(Area_q)))

## Now lets get the sample IVT value

sample_IVT = 1/9.8 * Area_q * Area_V / 800

print("\n\nSo, the IVT value for this vertical column assumed to be 1 x 1 m^2 is (1/g)(Area_q * Area_V)/dx to account for multiplying twice")
print("For this sample the IVT value is: {} kg/ms".format(sample_IVT))

class IVT:
    IVT = np.array([])
    def __init__(self,q,V,time):
        #Going to construct a numpy array of IVT values
        self.q_arr = q
        self.V_arr = V
        self.time = time
    def calcIVT_Trapz(self):
        dP=800
        #Calculate the IVT
        
        q_Area = np.array(np.abs([trapz(self.q_arr[t], dx = dP, axis = 0) for t in range(len(self.time))]))
        V_Area = np.array(np.abs([trapz(self.V_arr[t], dx = dP, axis = 0) for t in range(len(self.time))]))
        
        IVT_Trapz = 1/9.8 * q_Area * V_Area / dP
        return IVT_Trapz
    def calcIVT_Simps(self):
        dP = 800
        #Calculate the IVT
        q_Area = np.array(np.abs([simps(self.q_arr[t], dx = dP, axis = 0) for t in range(len(self.time))]))
        V_Area = np.array(np.abs([simps(self.V_arr[t], dx = dP, axis = 0) for t in range(len(self.time))]))
       
        IVT_Simps = 1/9.8 * q_Area * V_Area / dP
        return IVT_Simps
var = IVT(spec_humidity,u_wind,time)
IVT_arr_Trapz = var.calcIVT_Trapz() 
IVT_arr_Simps = var.calcIVT_Simps()
