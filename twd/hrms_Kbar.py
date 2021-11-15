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

sub_bx=int(sys.argv[9]) # Subdomain x index
sub_by=int(sys.argv[10]) # Subdomain y index

# --- SET PARAMETERS

# set FIGURE info
figdir  = workdir+'/plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

# set BOXES info for GEBCO 1/30'  
resol     = 120 # Npt/deg = 1'
boxdim    = 5   # deg

# set SUBDOMAINS info
[by,bx] = [sub_by,sub_bx]


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
print ('... mesh created. (NX,NY)=', NX,NY)

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

# Med subdivision in 48 subdomains
DOMY = [0,480,720,961,1200,NY-1] # 240
DOMX = [0,840,1680,2520,3360,4200,5040,5880,NX-1]

# MAIN CALC
print ('Start computation :')
print ('Compute running mean over 5x5deg boxes with steps of 1 grid point ...')
start = time.time()

for ji in range(DOMX[bx],DOMX[bx+1]) :
   #print ('Running index: ',ji)
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
   for jj in range(DOMY[by],DOMY[by+1]) :
      #print ('Running indexes: ',ji,jj)
      # all points near south and north boundaries will be equal to their interior neighbours
      [ya,yb] = [jj-mid, jj+mid] 

      # cut the box mask
      if centre_reg :
         #print ('CASE centre_reg msk[ya:yb,xa:xb]', msk[ya:yb,xa:xb])
         msk_b = msk[ya:yb,xa:xb]
      elif left_reg or right_reg :
         #print ('CASE left_reg or right_reg')
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
         #print ('hrms ',hrms_b)
         # calc K_bar
         integ  = np.sum( Kmod*hhat2*dk*dl ,axis=(0,1))
         frac   = 4*Area*(np.pi*hrms_b)**2
         Kbar_b = integ/frac
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

# --- write NPY file
temp_hrms_outfile = workdir+'/'+npy_hrms_pre+'dom'+str(by)+str(bx)
print('Saving: [%s]' %(temp_hrms_outfile+'.npy'))
np.save(temp_hrms_outfile,np.asarray(h_rms))
temp_kbar_outfile = workdir+'/'+npy_kbar_pre+'dom'+str(by)+str(bx)
print('Saving: [%s]' %(temp_kbar_outfile+'.npy'))
np.save(temp_kbar_outfile,np.asarray(K_bar))

