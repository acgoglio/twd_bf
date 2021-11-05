'''
Written by Federica Borile , Apr 2021
Modified by AC Goglio , Nov 2021
'''

print ("\n $$$--- CALC and SAVE h_rms and Kbar interpolated from GEBCO bathymetry to MED grid ---$$$ ")

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
bathy_inlat=str(sys.argv[6]) # lat field name in the input/output file
bathy_inlon=str(sys.argv[7]) # lon field name in the input/output file


# --- SET PARAMETERS

# set FIGURE info
figdir  = './plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

# set BOXES info for GEBCO 1/30'  
resol     = 120 # Npt/deg = 1'
boxdim    = 5   # deg

# set SUBDOMAINS info
# bx: 0,1,2
# by: 0,1,2,3
[by,bx] = [0,0]

# set FLAG
flag_hrms_plot = True
flag_Kbar_plot = True
flag_vars_save = True


#### ------- DO NOT CHANGE below this line ------- #####

def find_point_distance(lonA,latA,lonB,latB) :
   # coord are in deg xx.xxx°
   R = 6371229          # Earth radius [m]
   rad = np.pi / 180.0  # conversion from degree into radians
   # conversion to rad: work with matrices
   lon1 = lonA * rad     
   lat1 = latA * rad
   lon2 = lonB * rad
   lat2 = latB * rad
   dlon = lon2 - lon1
   dlat = lat2 - lat1
   a  = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
   c  = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
   dd = R * c   # distance in [m]
   return dd

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
ncf.close()

# Compute the  mesh
print ('Create the mesh ...')
[NY,NX]    = rough.shape
[nav_lon,nav_lat] = np.asarray(np.meshgrid(lon,lat)) 
e1t = find_point_distance(nav_lon[:,0:-1],nav_lat[:,0:-1],nav_lon[:,1:],nav_lat[:,1:])
e2t = find_point_distance(nav_lon[0:-1,:],nav_lat[0:-1,:],nav_lon[1:,:],nav_lat[1:,:])
# since it is an array of differences it has 1 point less respect to the grid
e1t = np.append(e1t,np.atleast_2d(e1t[:,-1]).T,axis=1)
e2t = np.append(e2t,np.atleast_2d(e2t[-1,:]),  axis=0)
area_cell = e1t * e2t
print ('... mesh created.', NY,NX)

#    CALC : as Shakespeare et al (2020) 
# ----------
# prepare variables
Npts  = boxdim*resol
mid   = Npts // 2
h_rms = np.empty(shape=(NY,NX))
K_bar = np.empty(shape=(NY,NX))
msk   = np.ones(shape=(NY,NX))
msk[bathy>=0.] = 0.
TM2 = 12.42 # hours, M2 tidal period
rad = np.pi / 180.0  # conversion from degree into radians
#latM2 =  np.arcsin(24/(2*TM2)) / rad
#jjmin = find_nearest(lat,-latM2) - 5
#jjmax = find_nearest(lat, latM2) + 5
#DOMX = [0    , 2700, 5400, 10800, 16200, NX-1]
#DOMY = [jjmin, 3900, 6900 , jjmax]
DOMX = [0,961,NX-1]
DOMY = [0,1680,3360,5040,NY-1]

# MAIN CALC
print ('Start computation :')
print ('Compute running mean over 5x5deg boxes with steps of 1 grid point ...')
start = time.time()

#for ji in range(DOMX[bx],DOMX[bx+1]) :
for ji in range(0,NX) :
   print ('Running index: ',ji)
   # define different x-position ranges to compose boxes near periodic boundaries
   [left_reg,centre_reg,right_reg] = [False, False, False]
   if ji<mid :
      left_reg   = True
      [xa,xb]    = [NX-(mid-ji)-1, ji+mid]
   elif ji>=mid and ji<(NX-1-mid) :
      centre_reg = True
      [xa,xb]    = [ji-mid, ji+mid]
   elif ji>=(NX-1-mid) :
      right_reg  = True
      [xa,xb]    = [ji-mid, ji-(NX-mid-1)]

   # compute only lats with f<M2freq (higher are not used in wave drag parametrization)
   #for jj in range(DOMY[by],DOMY[by+1]) :
   for jj in range(0,NY) :
      print ('Running indexes: ',ji,jj)
      # all points near south and north boundaries will be equal to their interior neighbours
      [ya,yb] = [jj-mid, jj+mid] 

      # cut the box mask
      if centre_reg :
         msk_b = msk[ya:yb,xa:xb]
      elif left_reg or right_reg :
         msk_b = np.concatenate((msk[ya:yb,xa:-1],msk[ya:yb,0:xb]), axis=1)
      else: 
         print ('ERROR!!!')
         sys.exit()
      # check ocean presence, box has AT LEAST ONE point 
      if np.nansum(msk_b)>0 :  
         # cut matricex boxes
         if centre_reg :
            Area = area_cell[ya:yb,xa:xb]
            h    = rough[ya:yb,xa:xb]
            dx   = e1t[ya:yb,xa:xb]
            dy   = e2t[ya:yb,xa:xb]
         elif left_reg or right_reg :
            Area = np.concatenate((area_cell[ya:yb,xa:-1],area_cell[ya:yb,0:xb]), axis=1)
            h    = np.concatenate((rough[ya:yb,xa:-1]    ,rough[ya:yb,0:xb])    , axis=1)
            dx   = np.concatenate((e1t[ya:yb,xa:-1]      ,e1t[ya:yb,0:xb])      , axis=1)
            dy   = np.concatenate((e2t[ya:yb,xa:-1]      ,e2t[ya:yb,0:xb])      , axis=1)
         Area = np.sum(Area,axis=(0,1)) # working with FT all points are considered even if they are dry!
         [Nyb,Nxb] = h.shape

         # calc mean dx and dy: we approximate them as const to compute FFT
         dx = np.average(dx, axis=(0,1))
         dy = np.average(dy, axis=(0,1))
         # calc FFT
         hhat  = np.fft.fft2(h)
         hhat  = hhat*dx*dy   # to have correct dim
         hhat2 = np.power(np.absolute(hhat),2)
         # calc (k,l) wavenumber
         fac = 2*np.pi
         k1d = np.fft.fftfreq(Nxb,d=dx) * fac
         l1d = np.fft.fftfreq(Nyb,d=dy) * fac
         dk  = np.abs(k1d[1]-k1d[0])
         dl  = np.abs(l1d[1]-l1d[0])
         [k,l] = np.asarray(np.meshgrid(k1d,l1d)) 
         Kmod  = np.sqrt(k**2+l**2)
         # calc h_rms
         integ  = np.sum( hhat2*dk*dl ,axis=(0,1)) # integral
         frac   = 4*Area*(np.pi)**2
         hrms_b = np.sqrt( integ/frac )
         # calc K_bar
         integ  = np.sum( Kmod*hhat2*dk*dl ,axis=(0,1))
         frac   = 4*Area*(np.pi*hrms_b)**2
         Kbar_b = integ/frac 
      else :
         hrms_b = 0.
         Kbar_b = 0.
      h_rms[jj,ji] = hrms_b
      K_bar[jj,ji] = Kbar_b
print ('... computation done.')
end = time.time()
print(end-start)

## apply mask to the final field (to be sure...)
#h_rms = np.where(msk==1, h_rms, -9999)
#K_bar = np.where(msk==1, K_bar, -9999)
## copy the first column on the last position
#h_rms[:,NX-1] = h_rms[:,0]
#K_bar[:,NX-1] = K_bar[:,0]

# --- write NPY file
outfile = 'hrms_temp_dom'+str(by)+str(bx)+'a' 
print('Saving: [%s]' %(outfile+'.npy'))
np.save(outfile,np.asarray(h_rms))
outfile = 'kbar_temp_dom'+str(by)+str(bx)+'a'
print('Saving: [%s]' %(outfile+'.npy'))
np.save(outfile,np.asarray(K_bar))


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
   latout = nco.createVariable('lat', dtype('double').char, 'y')
   latout.standard_name = 'latitude'
   latout.long_name = 'latitude'
   latout.units = 'degrees_north'
   latout.axis = 'Y'
   # create longitude axis
   lonout = nco.createVariable('lon', dtype('double').char, 'x')
   lonout.standard_name = 'longitude'
   lonout.long_name = 'longitude'
   lonout.units = 'degrees_east'
   lonout.axis = 'X'
   # copy axis from original dataset
   lonout[:] = lon[:]
   latout[:] = lat[:]
   # create variable array : var1
   out_var1 = nco.createVariable('h_rms', dtype('double').char, ('y','x'))
   out_var1.long_name   = 'root mean square topograpthic height (roughness)'
   out_var1.description = 'h_rms = \sqrt( \int |FFT(roughness)|^2 dk dl / (4*pi^2*A) )'
   out_var1.units       = 'm'
   # create variable array : var2
   out_var2 = nco.createVariable('K_bar', dtype('double').char, ('y','x'))
   out_var2.long_name   = 'height-weightened-mean wavenumber of topograpthic height (roughness)'
   out_var2.description = 'K_bar =  \int (k^2+l^2)*|FFT(roughness)|^2 dk dl / (h_rms^2*4*pi^2*A)'
   out_var2.units       = 'm-1'
   # copy axis from original dataset
   out_var1[:]   = h_rms
   out_var2[:]   = K_bar
   # create global attributes
   nco.description = 'h_rms and K_bar from roughness (r) on running boxes of '+str(boxdim)+'x'+str(boxdim)+" deg, where r is original-smoothed GEBCO 1/30' bathymetry (smoothing has been applied with 2 iterations of 2nd order Shapiro filter"
   nco.reference = 'Shakespeare et al., 2020 : The Drag on the Barotropic Tide due to the Generation of Baroclinic Motion (JPO)'
   # close files
   nco.close()
   print ("Saving: [%s]" %outfile)


#   -----------------------
#   |    make the plot    |
#   -----------------------

# Area
# Mask invalid values
nav_lon[nav_lon==0]=np.nan
#VAR[VAR==0]=np.nan
# Cut requested area
nav_lat[np.logical_or(nav_lat < -89,nav_lat > 89)] = np.nan
nav_lon[np.logical_or(nav_lon < -179,nav_lon > 179)] = np.nan
# Mask variables
nav_lat = np.ma.masked_invalid(nav_lat)
nav_lon = np.ma.masked_invalid(nav_lon)


   # --- PLOT

if flag_hrms_plot :
   VAR = h_rms
   VARunit = '[m]'   
   VAR = np.ma.masked_invalid(VAR)

   reg_n = 'Med Sea'
   print('\nRegion: ' + reg_n +'\n')
   figname = figdir +'map_hrms_runbox'+str(boxdim)+'_fft'
   figtitle = r'$h_{rms} = \sqrt{\frac{1}{4A\pi^2} \int\int |\hat{h}|^2 dk dl}$'
   figname = figname +'_'+ reg_n +'.png'
   cmap        = plt.cm.gist_heat_r # Colormap
   [cmin,cmax] = [0,600]            # color min and max values
      
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_bnd[0],urcrnrlat=lat_bnd[1],llcrnrlon=lon_bnd[0],urcrnrlon=lon_bnd[1],resolution='i')
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

   reg_n = 'Med Sea'
   print('\nRegion: ' + reg_n +'\n')
   figname = figdir +'map_Kbar_runbox'+str(boxdim)+'_fft'
   figtitle = r'$\overline{K} = \frac{1}{4A\pi^2h_{rms}^2} \int\int (k^2+l^2) |\hat{h}|^2 dk dl$'
   figname = figname +'_'+ reg_n +'.png'
   cmap        = plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [0.,1.e-4]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_bnd[0],urcrnrlat=lat_bnd[1],llcrnrlon=lon_bnd[0],urcrnrlon=lon_bnd[1],resolution='i')
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
