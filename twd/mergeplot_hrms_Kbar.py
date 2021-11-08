'''
Fede, May 2021
'''

print ("\n $$$--- SAVE h_rms and Kbar  ---$$$ ")
print ("\n $$$--- Merge subdomains together (h_rms and K_bar in the GEBCO grid at each point using subdomains) ... ")

###############
import datetime
import time
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from pylab import *
from matplotlib import pyplot as plt
#import cmocean
import numpy as np
import math
import sys
import os
###############
# --- READ LINE ARGS

workdir=str(sys.argv[1]) # Work directory
bathy_infile=str(sys.argv[2]) # Name of the bathymetry file (includes original bathy field and roughness field)
fn=workdir+'/'+bathy_infile

outfile=workdir+'/'+str(sys.argv[3]) # Output file name

bathy_inname=str(sys.argv[4]) # original bathymetry field name in the input file
bathy_rough=str(sys.argv[5]) # roughness field name in the input/output file
bathy_inlat=str(sys.argv[6]) # lat field name in the input file
bathy_inlon=str(sys.argv[7]) # lon field name in the input file

sub_bx_num=int(sys.argv[8]) # Num of Subdomains in x 
sub_by_num=int(sys.argv[9]) # Num of Subdomains in y

out_hrms_name=str(sys.argv[10]) # hrms field name in outfile
out_kbar_name=str(sys.argv[11]) # kbar field name in outfile
out_lat_name=bathy_inlat  # lat field name in outfile
out_lon_name=bathy_inlon  # lon field name in outfile

# --- SET PARAMETERS

# set FIGURE info
figdir  = './plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

# set BOXES info for GEBCO 1/30'  
resol     = 120 # Npt/deg = 1'
boxdim    = 5   # deg

# set FLAG
flag_hrms_plot = True
flag_Kbar_plot = True
flag_vars_save = True


#### ------- DO NOT CHANGE below this line ------- #####

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# Read GEBCO bathymetry file
print ('Reading: [%s]' %fn)
ncf = Dataset(fn,mode='r')
bathy = ncf.variables[bathy_inname][:]  ; bathy = np.squeeze(bathy)
lat   = ncf.variables[bathy_inlat][:]        ; lat = np.squeeze(lat)    # Y-axis
lon   = ncf.variables[bathy_inlon][:]        ; lon = np.squeeze(lon)    # X-axis
lat_bnd=[np.min(lat),np.max(lat)]
lon_bnd=[np.min(lon),np.max(lon)]
ncf.close()

# Read GEBCO roughness file
print ('Reading: [%s]' %fn)
ncf = Dataset(fn,mode='r')
rough = ncf.variables[bathy_rough][:] ; rough = np.squeeze(rough)
[NY,NX]    = rough.shape
ncf.close()

#    CALC : as Shakespeare et al (2020) 
# ----------
# prepare variables
hrms_f = np.empty(shape=(NY,NX))
kbar_f = np.empty(shape=(NY,NX))
msk   = np.ones(shape=(NY,NX))
msk[bathy>=0.] = 0.
TM2 = 12.42 # hours, M2 tidal period
rad = np.pi / 180.0  # conversion from degree into radians
latM2 =  np.arcsin(24/(2*TM2)) / rad
jjmin = find_nearest(lat,-latM2) - 5
jjmax = find_nearest(lat, latM2) + 5

# Med subdivision in 8 subdomains
DOMX = [0,961,NX-1]
DOMY = [0,1680,3360,5040,NY-1]

# Merge subdomains
for bx in range(0,sub_bx_num+1):
    for by in range(0,sub_by_num+1):
        subdom = str(bx)+str(by)
        print('I am merging the subdomain: ',subdom)

        tmp_h = np.load(workdir+'/'+npy_hrms_pre+'dom'+subdom+'.npy')
        tmp_k = np.load(workdir+'/'+npy_kbar_pre+'dom'+subdom+'.npy')
        hrms_f[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]] = tmp_h[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]] 
        kbar_f[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]] = tmp_k[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]]

# apply mask to the final field (to be sure...)
hrms_f = np.where(msk==1, hrms_f, -9999.)
kbar_f = np.where(msk==1, kbar_f, -9999.)


#   ---------------------
#   | write netCDF file |
#   --------------------

if flag_vars_save :
   # open netCDF file to write
   nco  = Dataset(outfile,mode='w',format='NETCDF4')

   # define axis size
   #ncout.createDimension('time', None)  # unlimited
   nco.createDimension('y', NY)
   nco.createDimension('x', NX)
   # create latitude axis
   latout = nco.createVariable(out_lat_name, dtype('double').char, 'y')
   latout.standard_name = 'latitude'
   latout.long_name = 'latitude'
   latout.units = 'degrees_north'
   latout.axis = 'Y'
   # create longitude axis
   lonout = nco.createVariable(out_lon_name, dtype('double').char, 'x')
   lonout.standard_name = 'longitude'
   lonout.long_name = 'longitude'
   lonout.units = 'degrees_east'
   lonout.axis = 'X'
   # copy axis from original dataset
   lonout[:] = lon[:]
   latout[:] = lat[:]
   # create variable array : var1
   out_var1 = nco.createVariable(out_hrms_name, dtype('double').char, ('y','x'))
   out_var1.long_name   = 'root mean square topograpthic height (roughness)'
   out_var1.description = 'h_rms = \sqrt( \int |FFT(roughness)|^2 dk dl / (4*pi^2*A) )'
   out_var1.units       = 'm'
   # create variable array : var2
   out_var2 = nco.createVariable(out_kbar_name, dtype('double').char, ('y','x'))
   out_var2.long_name   = 'height-weightened-mean wavenumber of topograpthic height (roughness)'
   out_var2.description = 'K_bar =  \int (k^2+l^2)*|FFT(roughness)|^2 dk dl / (h_rms^2*4*pi^2*A)'
   out_var2.units       = 'm-1'
   # copy axis from original dataset
   out_var1[:]   = hrms_f
   out_var2[:]   = kbar_f
   # create global attributes
   nco.description = 'h_rms and K_bar from roughness (r) on running boxes of '+str(boxdim)+'x'+str(boxdim)+" deg, where r is original-smoothed GEBCO 1' bathymetry (smoothing has been applied with 200 iterations of 2nd order Shapiro filter"
   nco.reference = 'Shakespeare et al., 2020 : The Drag on the Barotropic Tide due to the Generation of Baroclinic Motion (JPO)'
   # close files
   nco.close()
   print ("Saving: [%s]" %outfile)


#   -----------------------
#   |    make the plot    |
#   -----------------------

## Area
## Mask invalid values
#nav_lon[nav_lon==0]=np.nan
##VAR[VAR==0]=np.nan
## Mask variables
#nav_lat = np.ma.masked_invalid(nav_lat)
#nav_lon = np.ma.masked_invalid(nav_lon)


   # --- PLOT

if flag_hrms_plot :
   VAR = h_rms
   VARunit = '[m]'   
   VAR = np.ma.masked_invalid(VAR)

   k = 0 #for k,reg_n in enumerate(region):
   reg_n = 'global'
   print('\nRegion: ' + reg_n +'\n')
   figname = figdir +'map_hrms_gebco1_runbox'+str(boxdim)+'_step025_fft'
   figtitle = r'$h_{rms} = \sqrt{\frac{1}{4A\pi^2} \int\int |\hat{h}|^2 dk dl}$'
   figname = figname +'_'+ reg_n +'.png'
   cmap        = plt.cm.gist_heat_r # Colormap
   [cmin,cmax] = [0,600]            # color min and max values
      
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_bnd[k,0],urcrnrlat=lat_bnd[k,1],llcrnrlon=lon_bnd[k,0],urcrnrlon=lon_bnd[k,1],resolution='i')
   #      m.drawparallels(np.arange(-90,90,20),labels=[1,0,0,0], fontsize=12, linewidth=0.0)
   #      m.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1], fontsize=12, linewidth=0.0)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   pc  = plt.contour(x,y,bathy, levels=[1000], colors='dimgray')
   m.fillcontinents(color='0.8',lake_color='0.9')
   m.drawcoastlines(color='dimgray', linewidth=0.3)
   plt.title( figtitle, fontsize='18')
   cbar = m.colorbar(fig,'right', size='3%', pad='2%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
   
   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')

if flag_Kbar_plot :
   VAR = K_bar
   VARunit = r'[$m^{-1}]'
   VAR = np.ma.masked_invalid(VAR)

   k = 0 #for k,reg_n in enumerate(region):
   reg_n = 'global'
   print('\nRegion: ' + reg_n +'\n')
   figname = figdir +'map_Kbar_gebco1_runbox'+str(boxdim)+'_step025_fft'
   figtitle = r'$\overline{K} = \frac{1}{4A\pi^2h_{rms}^2} \int\int (k^2+l^2) |\hat{h}|^2 dk dl$'
   figname = figname +'_'+ reg_n +'.png'
   cmap        = plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [0.,1.e-4]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_bnd[k,0],urcrnrlat=lat_bnd[k,1],llcrnrlon=lon_bnd[k,0],urcrnrlon=lon_bnd[k,1],resolution='i')
   #      m.drawparallels(np.arange(-90,90,20),labels=[1,0,0,0], fontsize=12, linewidth=0.0)
   #      m.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1], fontsize=12, linewidth=0.0)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   pc  = plt.contour(x,y,bathy, levels=[1000], colors='dimgray')
   m.fillcontinents(color='0.8',lake_color='0.9')
   m.drawcoastlines(color='dimgray', linewidth=0.3)
   plt.title( figtitle, fontsize='18')
   cbar = m.colorbar(fig,'right', size='3%', pad='2%', extend='max') #, format='%.0e')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
   cbar.formatter.set_powerlimits((0, 0))

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')



#
## --- END
#
