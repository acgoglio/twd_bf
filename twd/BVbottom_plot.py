"""
Written by Federica Borile , May 2020
Modified by AC Goglio , Nov 2021
"""

print ('\n $$$--- Compute and SAVE the Brunt-Vaisala Frequency at the bottom in EAS Med system ---$$$\n')

###############
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from pylab import *
from matplotlib import pyplot as plt
import matplotlib.colors as colors
#import cmocean
import numpy as np
import sys
import os
###############
# --- READ LINE ARGS

workdir=str(sys.argv[1])   # Work directory
eas_bathy=str(sys.argv[2]) # Path/Name of the EAS system bathymetry file
eas_mesh=str(sys.argv[3])  # Path/Name of the EAS system mesh mask file

infile=workdir+'/'+str(sys.argv[4])    # Input file name
infield=str(sys.argv[5])   # Name of the mean BV field in the input file  

outfile=workdir+'/'+str(sys.argv[6])   # Output file name
outfield=str(sys.argv[7])   # Name of the bottom mean BV field in the ouput file  

eas_mbathy=str(sys.argv[8])        # eas mbathy field name in the mesh mask file
eas_tmask=str(sys.argv[9])         # eas tmask field name in the mesh mask file
eas_lat=str(sys.argv[10])           # eas lat field name in the mesh mask file
eas_lon=str(sys.argv[11])           # eas lon field name in the mesh mask file
eas_Bathymetry=str(sys.argv[12])    # eas Bathymetry field name in the bathy file

order=int(sys.argv[13]) # Order of the Shapiro filter
napp=int(sys.argv[14]) # Number of Shapiro filter applications
scheme=int(sys.argv[15]) # type of boundary scheme to use ( only option 1 implemented = No change at wall, constant order )  


# --- SET PARAMETERS

# set FIGURE info
figdir  = workdir+'/plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

# set FLAG
flag_calc_bnbot = True
flag_save_bnbot = True
flag_plot_shap   = True


##### ------- DO NOT CHANGE below this line ------- #####

# --- SET PARAMETERS
# Read mesh info
filemesh = eas_mesh
print ("Reading: [%s]" %filemesh)
mesh = Dataset(filemesh,mode='r')
nav_lat = mesh.variables[eas_lat][:]   ; nav_lat = np.squeeze(nav_lat)        # Y-axis
nav_lon = mesh.variables[eas_lon][:]   ; nav_lon = np.squeeze(nav_lon)        # X-axis
mbathy  = mesh.variables[eas_mbathy][:]    ; mbathy  = np.squeeze(mbathy)
tmask   = mesh.variables[eas_tmask][:]     ; tmask   = np.squeeze(tmask)
mesh.close()
[NZ,NY,NX] = tmask.shape

# Read bathymetry file
filebathy = eas_bathy
print ('Reading: [%s]' %filebathy)
ncf = Dataset(filebathy,mode='r')
bathy = ncf.variables[eas_Bathymetry][:]   ; bathy = np.squeeze(bathy)
ncf.close()


# --------- extract "bnb" from "bn"
if flag_calc_bnbot :
   # read mean bn from input file
   # note : bn is on w-points, so the corresponding mask is Wmask equalt to Tmask
   fn = infile
   print ("Reading: ",fn)
   fnc = Dataset(infile,mode='r')
   bn = fnc.variables[infield][:]  ; bn = np.squeeze(bn)
   #
   # extract bn at the bottom
   bn = bn*tmask
   bnb = np.empty([NY,NX])
   for i in range(NX) :
      for j in range(NY) :
         bmk = mbathy[j,i]-1
         if bmk>=0 :
            bnb[j,i] = bn[bmk,j,i]

# --------- define the Shapiro filter function
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

# --------- apply the filtering to bnbot field
# 2D Filtering of the field
F=bnb
Fout=F
Im=NX
Jm=NY
print ('Grid dims are: (lon,lat)= ',Im,Jm)
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
   bnb=F

# --------- save "bnb" in netCDF file
if flag_save_bnbot :
   var = bnb[:]
   var_name     = outfield
   var_longname = 'Buoyancy frequency at the bottom'
   var_details  = 'BVF computed at the bottom using data from 60 days of simulation without tides'
   var_units    = 's-1'
   # open a netCDF file to write
   outf = outfile
   nco  = Dataset(outf,mode='w', format='NETCDF4')
   # define axis size
   nco.createDimension('y', var.shape[0])
   nco.createDimension('x', var.shape[1])
   # create latitude axis
   lat = nco.createVariable(eas_lat, dtype('double').char, ('y','x'))
   lat.standard_name = 'latitude'
   lat.long_name = 'latitude'
   lat.units = 'degrees_north'
   lat.axis = 'Y'
   lat[:] = nav_lat[:]
   # create longitude axis
   lon = nco.createVariable(eas_lon, dtype('double').char, ('y','x'))
   lon.standard_name = 'longitude'
   lon.long_name = 'longitude'
   lon.units = 'degrees_east'
   lon.axis = 'X'
   lon[:] = nav_lon[:]
   # create variable array
   vout = nco.createVariable(var_name, dtype('double').char, ('y','x'))
   vout.long_name = var_longname
   vout.units     = var_units
   vout.details   = var_details
   vout[:] = var[:]
   # create mask array
   mout = nco.createVariable('tmask', dtype('double').char, ('y','x'))
   mout[:] = tmask[0,:,:]
   # close files
   nco.close()
   print ("Saving: [%s]" %outf)

# --------
# !!! now shapiro filter should be applied using MATLAB
# !!! so the filtered BNBOT should be read from a different file 
# --------


if flag_plot_shap :
   # Read file
   fn = outfile
   print ('Reading: [%s]' %fn)
   ncf = Dataset(fn,mode='r')
   bnbot = ncf.variables[outfield][:]   
   ncf.close()

   # Area
   # Mask invalid values
   nav_lon[nav_lon==0]=np.nan
   # Mask variables
   nav_lat = np.ma.masked_invalid(nav_lat)
   lat_min=np.min(nav_lat[:,0])
   lat_max=np.max(nav_lat[:,0])
   lon_min=np.min(nav_lon[0,:])
   lon_max=np.max(nav_lon[0,:])
   nav_lon = np.ma.masked_invalid(nav_lon)
   
   # --- PLOT
   VAR = bnbot
   VARunit = r'[$s^{-1}$]'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   k = 0 #for k,reg_n in enumerate(region):
   figname = figdir +'map_bottomBV.png'
   figtitle = r'$N_{bottom}$'
   cmap        = cm.get_cmap('bone_r') # Colormap
   cmap        = cm.get_cmap('twilight') # Colormap
   #[cmin,cmax] = [1.e-4,VAR.max()]            # color min and max values
   #[cmin,cmax] = [1.e-8,VAR.max()]            # color min and max values
   [cmin,cmax] = [1.e-5,1.e+0]
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   #if k!=0 :
   #   m.drawparallels(np.arange(-18,36,10),labels=[1,0,0,0], fontsize=12, linewidth=0.3)
   #   m.drawmeridians(np.arange(-18,36,10),labels=[0,0,0,1], fontsize=12, linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, norm=colors.LogNorm(vmin=cmin, vmax=cmax) )
   pc  = plt.contour(x,y,bathy, levels=[1000], colors='dimgray')
   m.fillcontinents(color='0.8',lake_color='0.9')
   m.drawcoastlines(color='dimgray', linewidth=0.3)
   plt.title( figtitle, fontsize='18')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
   
   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')


#
## --- END
#
