# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:19:55 2020

@author: willk
"""
import numpy as np
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

#Task 2: Build Data Structure###################################################################################################
    #Subtask: Build function to calculate IVT
#def IVT(q,V):# IVT = (1/g) * S_1000^300 q V dP
#    g = -9.8
#    (1/g) * (-700) * q * 700
#    #Main Task
u_wind_dict = {}
spec_humid_dict = {}
ivt_dict = {}
dict = {'u_wind':u_wind_dict, 'specific_humidity':spec_humid_dict, 'IVT': ivt_dict}

for i in range(len(time)):
    #u_wind_dict.update({time[i]:{}})
    for keys in dict:
        dict[keys].update({time[i]:{}})
for times in dict['u_wind']:
    for i in range(len(levs)): #Tried to do this in one line. Cannot hash iterable item
        u_wind_dict[times].update({levs[i]:np.array([u_wind[int(times)][i]])})
        spec_humid_dict[times].update({levs[i]:np.array([spec_humidity[int(times)][i]])})
        #ivt_dict[times].update({levs[i]:np.array([])}) Will need to calculate this later

     
    

#Task 3: Build Map ###################################################################################################
lon_0 = np.mean(lons)
lat_0 = np.mean(lats)    
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')


#Task 4: plotv###################################################################################################
for times in dict['u_wind']:
    #print(yearmonth)
    for levels in dict['u_wind'][times]:
        #fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        #if (heights_dict[yearmonth][days] > 0.0) > 0.0:
        #    data = np.array(heights_dict[yearmonth][days])
        #data = np.array(heights_dict[yearmonth][days])
        #print(yearmonth,days)
        
        #    print("true")
            #N_heights_dict[yearmonth][days] = N_heights_dict[yearmonth][days].filled(np.mean(N_heights_dict[yearmonth][days]))
            #cs = m.pcolor(xi,yi,np.squeeze(N_heights_dict[yearmonth][days]), cmap = 'hsv')
        #u_wind
        cs1 = m.pcolor(lons, lats , np.squeeze(u_wind_dict[times][levels]), cmap = 'hsv')
        m.drawparallels(np.arange(-90.,91.,30.))
        m.drawmeridians(np.arange(-180.,181.,60.))

        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()


        cbar = m.colorbar(cs1, location='right', pad="10%")
        cbar.set_label('m/s')
                
        ax1.set_title('u_wind at {} hpa: Day: {}'.format(levels, times))
        
        
        #Specific Humidity
        ax2 = fig.add_subplot(222)
        cs2 = m.pcolor(lons, lats , np.squeeze(spec_humid_dict[times][levels]), cmap = 'RdYlGn')
        m.drawparallels(np.arange(-90.,91.,30.))
        m.drawmeridians(np.arange(-180.,181.,60.))

        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()


        cbar = m.colorbar(cs2, location='right', pad="10%")
        cbar.set_label('g/kg')
        
        ax2.set_title('Specific humidity at {} hpa: Day: {}'.format(levels, times))
        plt.tight_layout()
        plt.show()
        #
        #ax3 = fig.add_subplot(223)
        #cs3 = m.pcolor(lons, lats , np.squeeze(ivt_dict[times][levels]), cmap = 'Blues')
        #m.drawparallels(np.arange(-90.,91.,30.))
        #m.drawmeridians(np.arange(-180.,181.,60.))

        #m.drawcoastlines()
        #m.drawstates()
        #m.drawcountries()

        #I have cahnge somehting here
        
        #cbar = m.colorbar(cs3, location='right', pad="10%")
        #cbar.set_label('kg/ms')
        
        #ax3.set_title('Integrated Vapor Transport at\n {} hpa: Day: {}'.format(levels, times))
        plt.tight_layout()
        plt.show()
    



#Task 5: Sorting Algorithm

 
