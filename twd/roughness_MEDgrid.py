# Import libraries
# imports
import matplotlib.pyplot as plt
import matplotlib as mpl # Palettes
import numpy as np
import netCDF4 as NC
import os
import sys
import warnings
#
warnings.filterwarnings("ignore") # Avoid warnings
#
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from datetime import date, timedelta
from datetime import datetime
from operator import itemgetter
import plotly
from plotly import graph_objects as go # for bar plot
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from operator import itemgetter # to order lists
from statsmodels.distributions.empirical_distribution import ECDF # empirical distribution functions

# Read the line args
workdir=str(sys.argv[1]) # Work directory
bathy_infile=str(sys.argv[2]) # Name of the input/output bathymetry file
bathy_inname=str(sys.argv[3]) # original bathymetry field name in the input/output file
bathy_outname=str(sys.argv[4]) # smoothed bathymetry field name in the input/output file
bathy_rough=str(sys.argv[5]) # roughness field name in the output file
bathy_inlat=str(sys.argv[6]) # lat field name in the input/output file
bathy_inlon=str(sys.argv[7]) # lon field name in the input/output file

eas_bathy=str(sys.argv[8]) # Path/Name of the EAS system bathymetry file
eas_Bathymetry=str(sys.argv[9])    # eas Bathymetry field name in the bathy file
eas_lat=str(sys.argv[10])           # eas lat field name in the mesh mask file
eas_lon=str(sys.argv[11])           # eas lon field name in the mesh mask file

np2be0=int(sys.argv[12]) # Num of Atlantic bdy points to be set to 0 (because of the artificialy modified EAS bathymetry) 

print ('You are running the script with the following line args: ', workdir,bathy_infile,bathy_inname,bathy_outname,bathy_rough,bathy_inlat,bathy_inlon,eas_bathy,eas_Bathymetry,eas_lat,eas_lon)

# Flags
flag_computesave_roughness=1
flag_plot_roughness=1

###########################
if flag_computesave_roughness == 1:

   # Open the input/output netCDF and read the original bathy field dims and values
   nc2open=workdir+'/'+bathy_infile 
   
   print ('I am going to open and modify the copy of the following file: ',nc2open)
   bathy_infield = NC.Dataset(nc2open,'r+')
   
   inbathymetry=bathy_infield.variables[bathy_inname][:]
   nav_lat = bathy_infield.variables[eas_lat][:]   ; nav_lat = np.squeeze(nav_lat)        # Y-axis
   nav_lon = bathy_infield.variables[eas_lon][:]   ; nav_lon = np.squeeze(nav_lon)        # X-axis
   
   # Create the new fields to be appended in the input/output netCDF
   print ('I am going to add the new bathymetry field: ',bathy_outname)
   outbathymetry=bathy_infield.createVariable(bathy_outname, np.float64, (bathy_inlat,bathy_inlon))
   outbathymetry.units = 'm'
   
   # Read the EAS bathymetry field
   print ('Reading: [%s]' %eas_bathy)
   bncf = NC.Dataset(eas_bathy,mode='r')
   easbathy = bncf.variables[eas_Bathymetry][0,:,:] 
   bncf.close()
   
   # Write the new fields in the input/output netCDF
   print ('I am writing the new field in the input/output netCDF..')
   print('easbathy',easbathy.shape)
   outbathymetry[:]=easbathy[:]
   
   bathy_infield.close()
   print ('Done')
   
   #######################
   # Add the differences
   # Open the input/output netCDF and read the original bathy field dims and values
   nc2open=workdir+'/'+bathy_infile
   
   print ('I am going to open and modify the following file: ',nc2open)
   bathy_infield = NC.Dataset(nc2open,'r+')
   
   # Read the fields
   inbathymetry=bathy_infield.variables[bathy_inname][:]
   outbathymetry=bathy_infield.variables[bathy_outname][:]
   
   # Create the new fields to be appended in the input/output netCDF
   print ('I am going to add the new bathymetry field: ',bathy_outname)
   diff_bathy=bathy_infield.createVariable(bathy_rough,np.float64,(bathy_inlat,bathy_inlon))
   diff_bathy.units = 'm'
   
   print ('I am going to compute and add the roughness field: ',bathy_rough)
   
   diff=inbathymetry[:]-outbathymetry[:]
   print ('Max diff', np.min(diff),np.max(diff),np.mean(diff))

   # Correction for Atlantic bdy: put roughness to 0 in the first np2be0 grid points
   # Set to 0 the first 10 points on the western bdy
   diff[0,:,0:np2be0] = 0.
   # Set to 0 the first 10 points on the lower bdy
   zero_lat_idx=int(np.where(nav_lon[0,:] == 0.)[0]+1)
   for x_idx in range (0,zero_lat_idx): 
       for y_idx in range (0,np2be0):
           if outbathymetry[y_idx,x_idx] != 0.:
              diff[0,y_idx,x_idx] = 0.
   # Set to 0 the first 10 points on the upper bdy
   NY=len(nav_lat[:,0])
   for x_idx in range (0,zero_lat_idx):
       for y_idx in range (NY-np2be0,NY):
           if outbathymetry[y_idx,x_idx] != 0.:
              diff[0,y_idx,x_idx] = 0.

   # Write the new field 
   diff_bathy[:]=diff[:]
   
   bathy_infield.close()
   print ('Done')

###############################
# Plot the Roughness field
if flag_plot_roughness == 1:

   # Input file path/name
   nc2open=workdir+'/'+bathy_infile 

   print ('I am going to open and plot the following file: ',nc2open)

   bathy_infield = NC.Dataset(nc2open,'r')

   nav_lat = bathy_infield.variables[eas_lat][:]   ; nav_lat = np.squeeze(nav_lat)        # Y-axis
   nav_lon = bathy_infield.variables[eas_lon][:]   ; nav_lon = np.squeeze(nav_lon)        # X-axis
   roughness=bathy_infield.variables[bathy_rough][:]
   inbathymetry=bathy_infield.variables[bathy_inname][0,:,:]

   bathy_infield.close()

   # Area
   # Mask variables
   nav_lat = np.ma.masked_invalid(nav_lat)
   lat_min=np.min(nav_lat[:,0])
   lat_max=np.max(nav_lat[:,0])
   lon_min=np.min(nav_lon[0,:])
   lon_max=np.max(nav_lon[0,:])
   nav_lon = np.ma.masked_invalid(nav_lon)
   print ('prova',lat_min,lat_max,lon_min,lon_max)

   # PLOT Abs values
   VAR = roughness
   LAND=inbathymetry
   VARunit = 'm'
   VAR = np.ma.masked_invalid(VAR)

   figdir  = workdir+'/plots/'
   if not(os.path.isdir(figdir)) :
      print('Creating: [%s]' %figdir)
      os.makedirs(figdir)

   figname = figdir +'map_med_rough.png'
   figtitle = 'Roughness'
   cmap        = 'seismic' #plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [-50,50]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rc('font', size=10)
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   #fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   fig = m.contourf(x,y,VAR, levels=[-90,-70,-50,-30,-10,10,30,50,70,90] ,cmap=cmap, extend='both') # levels=[-90,-70,-50,-30,-10,10,30,50,70,90]
   pcf  = plt.contourf(x,y,LAND, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,LAND, levels=[15.0], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='8')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label('Roughness ['+VARunit+']',fontsize='8') 
   cbar.ax.tick_params(labelsize='8')
   cbar.formatter.set_powerlimits((0, 0))

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')




   # PLOT rel values
   VAR = roughness/inbathymetry
   VAR = VAR * 100.0
   LAND=inbathymetry
   VARunit = ' ' 
   VAR = np.ma.masked_invalid(VAR)

   figdir  = workdir+'/plots/'
   if not(os.path.isdir(figdir)) :
      print('Creating: [%s]' %figdir)
      os.makedirs(figdir)

   figname = figdir +'map_med_rough_rel.png'
   figtitle = 'Roughness/Bathymetry'
   cmap        = 'seismic' #plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [-10,10]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rc('font', size=8)
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   #fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   fig = m.contourf(x,y,VAR, levels=[-9.0,-7.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,7.0,9.0] ,cmap=cmap, extend='both') # levels=[-90,-70,-50,-30,-10,10,30,50,70,90]
   pcf  = plt.contourf(x,y,LAND, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,LAND, levels=[15.0], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='8')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label('Roughness/Bathymetry [%]',fontsize='8') 
   cbar.ax.tick_params(labelsize='8')
   cbar.formatter.set_powerlimits((0, 0))

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')



