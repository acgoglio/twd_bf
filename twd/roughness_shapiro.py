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

order=int(sys.argv[8]) # Order of the Shapiro filter
napp=int(sys.argv[9]) # Number of Shapiro filter applications
scheme=int(sys.argv[10]) # type of boundary scheme to use ( only option 1 implemented = No change at wall, constant order )  

print ('You are running the script with the following line args: ', workdir,bathy_infile,bathy_inname,bathy_outname,bathy_inlat,bathy_inlon,order,napp,scheme)

# Open the input/output netCDF and read the original bathy field dims and values
nc2open=workdir+'/'+bathy_infile 

print ('I am going to open and modify the following file: ',nc2open)
bathy_infield = NC.Dataset(nc2open,'r+')

inbathymetry=bathy_infield.variables[bathy_inname][:]
lat=bathy_infield.variables[bathy_inlat][:]
lon=bathy_infield.variables[bathy_inlon][:]

Im=len(lon)
Jm=len(lat)
print ('Grid dims are: (lon,lat)= ',Im,Jm)

F=inbathymetry

# Create the new fields to be appended in the input/output netCDF
print ('I am going to add the new bathymetry field: ',bathy_outname)
outbathymetry=bathy_infield.createVariable(bathy_outname,np.float64,(bathy_inlat,bathy_inlon))
outbathymetry.units = 'm'

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
Fout=F
for n in range (1,napp+1):
   print(n,'^ application of the Shapiro filter ongoing..')

   # ----------------------------------------------------------------------------
   #  Filter all rows.
   # ----------------------------------------------------------------------------
   print ('I am going to filter the rows..')
   for j in range (0,Im):
       #print ('Filtering row num: ',j)
       Fraw=np.squeeze(F[:,j])
       # Run Shapiro 1D
       Fwrk=shapiro1D(Fraw,order,scheme)
       Fout[:,j]=Fwrk
       #print ('row Done!')

   # ----------------------------------------------------------------------------
   #  Filter all columns.
   # ----------------------------------------------------------------------------
   print ('I am going to filter the columns..')
   for i in range (0,Jm):
       #print ('Filtering col num: ',i)
       Fraw=np.squeeze(Fout[i,:])
       # Run Shapiro 1D
       Fwrk=shapiro1D(Fraw,order,scheme)
       Fout[i,:]=Fwrk
       #print ('row Done!')

   F=Fout

# Write the new fields in the input/output netCDF
print ('I am writing the new field in the input/output netCDF..')
outbathymetry[:]=np.squeeze(F[:])

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
lat=bathy_infield.variables[bathy_inlat][:]
lon=bathy_infield.variables[bathy_inlon][:]

# Create the new fields to be appended in the input/output netCDF
print ('I am going to add the new bathymetry field: ',bathy_outname)
diff_bathy=bathy_infield.createVariable(bathy_rough,np.float64,(bathy_inlat,bathy_inlon))
diff_bathy.units = 'm'

print ('I am going to compute and add the roughness field: ',bathy_rough)

diff=inbathymetry[:]-outbathymetry[:]
print ('Max diff', np.min(diff),np.max(diff),np.mean(diff))
diff_bathy[:]=diff[:]

bathy_infield.close()
print ('Done')

