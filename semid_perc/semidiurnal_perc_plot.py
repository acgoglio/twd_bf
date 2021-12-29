# imports
from __future__ import print_function, division
import os
import time
import sys
import warnings
#
warnings.filterwarnings("ignore") # Avoid warnings
#
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from operator import itemgetter
#
from numpy import ma
from scipy.stats import ks_2samp
from statsmodels.distributions.empirical_distribution import ECDF
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from scipy import interpolate
#
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 29/12/2021
# Modified: 29/12/2021
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# READ General run parameters: 
#---------------------
# Work directory
workdir_path=str(sys.argv[1])
# Name of the input file 
infile=str(sys.argv[2])
# Name of the field to be plotted
infield=str(sys.argv[3])
# Lon field name
inlon=str(sys.argv[4])
# Lat field name
inlat=str(sys.argv[5])

# Bathymetry file name and field name
model_bathy=str(sys.argv[6])
field_bathy=str(sys.argv[7])

# Mesh mask ad tmask field name
model_meshmask=str(sys.argv[8])
field_mesh=str(sys.argv[9])

# Name of the plot
outplot='SD_perc.png'
outplot=workdir_path+'/'+outplot

# ----------------------

# set FIGURE info
figdir  = workdir_path+'/plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)


# Build the path/name of the nc file and open it 
nc2open=workdir_path+infile
print ('Input file = ',nc2open)
model = NC.Dataset(nc2open,'r')
ST_perc=model.variables[infield][:]
vals=ST_perc

# Read lat, lon and fild values 
nc2open3=model_bathy # tidal bathimetry
model3 = NC.Dataset(nc2open3,'r')
nc2open4=model_meshmask # mesh mask
model4 = NC.Dataset(nc2open4,'r')

vals_bathy=model3.variables[field_bathy][:]
lons = model3.variables[inlon][:]
lats = model3.variables[inlat][:]

vals_land=model4.variables[field_mesh][0,0,:,:]
vals_land=np.squeeze(vals_land)

vals_ma = np.ma.masked_where(vals_land, vals)

plt.figure(figsize=(20,10))
plt.rc('font', size=16)

plt.title ('Semidiurnal Component Amplitude Percentages')

# Read the coordinates for the plot 
lon_0 = lons.mean()
llcrnrlon = lons.min()
urcrnrlon = lons.max()
lat_0 = lats.mean()
llcrnrlat = lats.min()
urcrnrlat = lats.max()

# Create the map
m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lons, lats)

# Plot the frame to the map
plt.rcParams["axes.linewidth"]  = 1.25

# Amp plot:
cs = plt.contourf(xi,yi,vals,levels=[0,10,20,30,40,50,60,70,80,90,100],cmap='coolwarm')
contour50 = plt.contour(xi,yi,vals,[50],colors='blue')

# Add the grid
m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)

# Land from bathymetry or mesh mask file
contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000],colors='black')
#contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,500.0], colors='gray')
#contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),[500.0],colors='black')

# Plot the legend 
cbar = m.colorbar(cs, location='bottom', pad="10%")
bar_label_string=' Semidiurnal component [%]'
cbar.set_label(bar_label_string)

plt.savefig(figdir+'/'+outplot)
print ('Outfile ',outplot)
plt.clf()

