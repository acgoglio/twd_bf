'''
Written by Federica Borile , Apr 2021
Modified by AC Goglio , Nov 2021
'''

print ("\n $$$--- CALC h_rms and Kbar interpolated from GEBCO bathymetry to MED grid ---$$$ ")

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

npy_hrms_pre=str(sys.argv[3]) # npy hrms Output file prename
npy_kbar_pre=str(sys.argv[4]) # npy kbar Output file prename

bathy_inname=str(sys.argv[5]) # original bathymetry field name in the input file
bathy_rough=str(sys.argv[6]) # roughness field name in the input/output file
bathy_inlat=str(sys.argv[7]) # lat field name in the input/output file
bathy_inlon=str(sys.argv[8]) # lon field name in the input/output file

out_hrms_name=str(sys.argv[9]) # hrms field name in both outfile
out_kbar_name=str(sys.argv[10]) # kbar field name in both outfile

outfile=str(sys.argv[11]) # Outfile name
outfile=workdir+'/'+outfile # Outfield name

resol = int(sys.argv[12]) # Resolution of the input fields
boxdim = int(sys.argv[13]) # Dimension of the box for the computation 

order_hk=int(sys.argv[14]) # Order of the Shapiro filter
napp_hk=int(sys.argv[15]) # Number of Shapiro filter applications
scheme_hk=int(sys.argv[16]) # type of boundary scheme to use ( only option 1 implemented = No change at wall, constant order ) 

# --- SET PARAMETERS

flag_eas_vars_save = True
flag_outfield_plot = True

# set FIGURE info
figdir  = workdir+'/plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

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
bathy = ncf.variables[bathy_inname][:]   ; bathy = np.squeeze(bathy)
nav_lat = ncf.variables[bathy_inlat][:]        # Y-axis
lat = nav_lat[:,0] 
lat = np.squeeze(lat)
nav_lon = ncf.variables[bathy_inlon][:]        # X-axis
lon = nav_lon[0,:]
lon = np.squeeze(lon)
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
print ('... mesh created. (NX,NY)=', NX,NY,'area cell= ',area_cell)

#    CALC : as Shakespeare et al (2020) 
# ----------
# prepare variables
Npts  = boxdim*resol
mid   = Npts // 2
h_rms = np.empty(shape=(NY,NX))
K_bar = np.empty(shape=(NY,NX))
msk   = np.ones(shape=(NY,NX))
msk[bathy<=0.] = 0.
TM2 = 12.42 # hours, M2 tidal period
rad = np.pi / 180.0  # conversion from degree into radians

# MAIN CALC
print ('Start computation :')
print ('Compute running mean over '+str(boxdim)+'x'+str(boxdim)+' deg boxes with steps of 1 grid point ...')
start = time.time()

for ji in range(0,NX):
   #print ('Running index: ',ji)
   # define different x-position ranges to compose boxes near boundaries
   # Left bdy:
   if ji<mid :
      [xa,xb]    = [0, ji+mid+(mid-ji)]
   # Internal 
   elif ji>=mid and ji<(NX-1-mid) :
      [xa,xb]    = [ji-mid, ji+mid]
   # Right bdy:
   elif ji>=(NX-1-mid) :
      [xa,xb]    = [ji-mid+(ji-(NX-mid-1)), NX-1]

   # compute only lats with f<M2freq (higher are not used in wave drag parametrization)
   for jj in range(0,NY) :
      #print ('Running indexes: ',ji,jj)
      # define different y-position ranges to compose boxes near boundaries
      # Lower bdy:
      if jj<mid :
         [ya,yb] = [0, jj+mid+(mid-jj)]
      # Internal
      elif jj>=mid and jj<(NY-1-mid) : 
         [ya,yb] = [jj-mid, jj+mid] 
      # Upper bdy:
      elif jj>=(NY-1-mid) :
         [ya,yb] = [jj-mid-(jj+mid-(NY-1)), NY-1]

      # cut the box mask
      msk_b = msk[ya:yb,xa:xb]

      # check ocean presence, box has AT LEAST ONE point 
      if np.nansum(msk_b)>0 :  
         # cut matricex boxes
         Area = area_cell[ya:yb,xa:xb]
         h    = rough[ya:yb,xa:xb]
         dx   = e1t[ya:yb,xa:xb]
         dy   = e2t[ya:yb,xa:xb]
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
         #print ('integ and frac',integ,frac)
         hrms_b = np.sqrt( integ/frac )
         #print ('hrms ',hrms_b)
         # calc K_bar
         if hrms_b != 0. :
            integ  = np.sum( Kmod*hhat2*dk*dl ,axis=(0,1))
            frac   = 4*Area*(np.pi*hrms_b)**2
            #print ('integ and frac',integ,frac)
            Kbar_b = integ/frac
         else:
            Kbar_b = 0.
         #print ('Kbar_b ',Kbar_b) 
      else :
         #print ('CASE: np.nansum(msk_b)<=0')
         hrms_b = 0.
         Kbar_b = 0.
      h_rms[jj,ji] = hrms_b
      K_bar[jj,ji] = Kbar_b
print ('... computation done.')
end = time.time()
print(end-start)

# -----------------------------------------------
# Define the Shapiro filter and smooth the fields

# 1D Shapiro filter function
# SHAPIRO 1D Function
def shapiro1D(Finp,order,scheme):

        Finp=np.array(Finp)
        order=np.array(order)
        scheme=np.array(scheme)

        ###########
        fourk=[2.500000E-1,6.250000E-2,1.562500E-2,3.906250E-3,9.765625E-4,2.44140625E-4,6.103515625E-5,1.5258789063E-5,3.814697E-6,9.536743E-7,2.384186E-7,5.960464E-8,1.490116E-8,3.725290E-9,9.313226E-10,2.328306E-10,5.820766E-11,1.455192E-11,3.637979E-12,9.094947E-13]

        Im1D=len(Finp)
        order2=int(np.fix(order/2))

        cor=[0 for i in range(Im1D)]
        Fcor=[0 for i in range(Im1D)]

        #----------------------------------------------------------------------------
        # Compute filter correction.
        #----------------------------------------------------------------------------

        if (scheme == 1):
           # Scheme 1:  constant order and no change at wall.

           # Filter loop
          for n in range (1,order2+1):

            # Set the bdy
            if (n != order2):
              cor[0]=2.0*(Finp[0]-Finp[1])
              cor[Im1D-1]=2.0*(Finp[Im1D-1]-Finp[Im1D-2])
            else:
              cor[0]=0.0
              cor[Im1D-1]=0.0

            # Set all the other
            cor[1:Im1D-2]=2.0*Finp[1:Im1D-2] - Finp[0:Im1D-3] - Finp[2:Im1D-1]

          coeff_to_apply=float(fourk[order2-1])
          #print ('Shapiro coeff. ',coeff_to_apply)
          Fcor=np.array(cor)*coeff_to_apply

        else:
          print ('Not yet implemented..')

        #----------------------------------------------------------------------------
        # Apply correction.
        #----------------------------------------------------------------------------

        Fout=Finp-Fcor

        return Fout

# 2D Filtering of the field
# hrms
Fout=h_rms
for n in range (1,napp_hk+1):
      print(n,'^ application of the Shapiro filter ongoing..')

      # ----------------------------------------------------------------------------
      #  Filter all rows.
      # ----------------------------------------------------------------------------
      print ('I am going to filter the rows..')
      for j in range (0,NX):
          #print ('Filtering row num: ',j)
          Fraw=np.squeeze(h_rms[:,j])
          # Run Shapiro 1D
          Fwrk=shapiro1D(Fraw,order_hk,scheme_hk)
          Fout[:,j]=Fwrk
          #print ('row Done!')

      # ----------------------------------------------------------------------------
      #  Filter all columns.
      # ----------------------------------------------------------------------------
      print ('I am going to filter the columns..')
      for i in range (0,NY):
          #print ('Filtering col num: ',i)
          Fraw=np.squeeze(Fout[i,:])
          # Run Shapiro 1D
          Fwrk=shapiro1D(Fraw,order_hk,scheme_hk)
          Fout[i,:]=Fwrk
          #print ('row Done!')

      h_rms=Fout

# Kbar
Fout=K_bar
for n in range (1,napp_hk+1):
      print(n,'^ application of the Shapiro filter ongoing..')

      # ----------------------------------------------------------------------------
      #  Filter all rows.
      # ----------------------------------------------------------------------------
      print ('I am going to filter the rows..')
      for j in range (0,NX):
          #print ('Filtering row num: ',j)
          Fraw=np.squeeze(K_bar[:,j])
          # Run Shapiro 1D
          Fwrk=shapiro1D(Fraw,order_hk,scheme_hk)
          Fout[:,j]=Fwrk
          #print ('row Done!')

      # ----------------------------------------------------------------------------
      #  Filter all columns.
      # ----------------------------------------------------------------------------
      print ('I am going to filter the columns..')
      for i in range (0,NY):
          #print ('Filtering col num: ',i)
          Fraw=np.squeeze(Fout[i,:])
          # Run Shapiro 1D
          Fwrk=shapiro1D(Fraw,order_hk,scheme_hk)
          Fout[i,:]=Fwrk
          #print ('row Done!')

      K_bar=Fout


# ------------------------------------
# Add the outfields in the outfile_eas
if flag_eas_vars_save :

   nc2open=outfile
   print ('I am going to open and modify the following file: ',nc2open)
   temp_file = Dataset(nc2open,'r+')

   # Create the new fields to be appended in the input/output netCDF
   print ('I am going to add the new hrms and kbar fields to the outfile: ',nc2open)
   new_eas_hrms=temp_file.createVariable(out_hrms_name,dtype('double').char, ('y','x'))
   new_eas_hrms.units = 'm'
   new_eas_kbar=temp_file.createVariable(out_kbar_name,dtype('double').char, ('y','x'))
   new_eas_kbar.units = 'm-1'

   # Add the fields
   new_eas_hrms[:]=h_rms
   new_eas_kbar[:]=K_bar

   temp_file.close()
   print ('Done')

# Plot the output fields
if flag_outfield_plot :

   # Input file path/name
   nc2open=outfile

   print ('I am going to open and plot the following file: ',nc2open)

   bathy_infield = Dataset(nc2open,'r')

   nav_lat = bathy_infield.variables[bathy_inlat][:]   ; nav_lat = np.squeeze(nav_lat)        # Y-axis
   nav_lon = bathy_infield.variables[bathy_inlon][:]   ; nav_lon = np.squeeze(nav_lon)        # X-axis
   roughness=bathy_infield.variables[bathy_rough][:]
   hrms_out=bathy_infield.variables[out_hrms_name][:]
   kbar_out=bathy_infield.variables[out_kbar_name][:]
   inbathymetry=bathy_infield.variables[bathy_inname][:]

   bathy_infield.close()

   # Area
   # Mask variables
   nav_lat = np.ma.masked_invalid(nav_lat)
   lat_min=np.min(nav_lat[:,0])
   lat_max=np.max(nav_lat[:,0])
   lon_min=np.min(nav_lon[0,:])
   lon_max=np.max(nav_lon[0,:])
   nav_lon = np.ma.masked_invalid(nav_lon)

   # PLOT HRMS
   VAR = hrms_out
   LAND=inbathymetry
   VARunit = 'm'
   VAR = np.ma.masked_invalid(VAR)

   figdir  = workdir+'/plots/'
   if not(os.path.isdir(figdir)) :
      print('Creating: [%s]' %figdir)
      os.makedirs(figdir)

   figname = figdir +'map_med_hrms_'+str(boxdim)+'x'+str(boxdim)+'_'+str(napp_hk)+'.png'
   figtitle = 'HRMS'
   cmap        = 'viridis' #plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [0,80]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rc('font', size=10)
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   #fig = m.contourf(x,y,VAR, levels=[-90,-70,-50,-30,-10,10,30,50,70,90] ,cmap=cmap, extend='both') # levels=[-90,-70,-50,-30,-10,10,30,50,70,90]
   pcf  = plt.contourf(x,y,LAND, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,LAND, levels=[15.0], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='8')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label('hrms ['+VARunit+']',fontsize='8')
   cbar.ax.tick_params(labelsize='8')
   cbar.formatter.set_powerlimits((0, 0))

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')

   # PLOT KBAR
   VAR = kbar_out
   LAND=inbathymetry
   VARunit = 'm-1'
   VAR = np.ma.masked_invalid(VAR)

   figdir  = workdir+'/plots/'
   if not(os.path.isdir(figdir)) :
      print('Creating: [%s]' %figdir)
      os.makedirs(figdir)

   figname = figdir +'map_med_kbar_'+str(boxdim)+'x'+str(boxdim)+'_'+str(napp_hk)+'.png'
   figtitle = 'Kbar'
   cmap        = 'magma' #plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [0,0.0008]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rc('font', size=8)
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   #fig = m.contourf(x,y,VAR, levels=[-9.0,-7.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,7.0,9.0] ,cmap=cmap, extend='both') # levels=[-90,-70,-50,-30,-10,10,30,50,70,90]
   pcf  = plt.contourf(x,y,LAND, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,LAND, levels=[15.0], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='8')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label('Kbar [m-1]',fontsize='8')
   cbar.ax.tick_params(labelsize='8')
   cbar.formatter.set_powerlimits((0, 0))

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')


