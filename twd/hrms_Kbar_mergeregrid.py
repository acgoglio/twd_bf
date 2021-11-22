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
from scipy import interpolate
###############
# --- READ LINE ARGS

workdir=str(sys.argv[1]) # Work directory
bathy_infile=str(sys.argv[2]) # Name of the bathymetry file (includes original bathy field and roughness field)
fn=workdir+'/'+bathy_infile

outfile_gebco=workdir+'/'+str(sys.argv[3]) # GEBCO-grid Output file name

bathy_inname=str(sys.argv[4]) # original bathymetry field name in the input file
bathy_rough=str(sys.argv[5]) # roughness field name in the input/output file
bathy_inlat=str(sys.argv[6]) # lat field name in the input file
bathy_inlon=str(sys.argv[7]) # lon field name in the input file

sub_bx_num=int(sys.argv[8]) # Num of Subdomains in x 
sub_by_num=int(sys.argv[9]) # Num of Subdomains in y

out_hrms_name=str(sys.argv[10]) # hrms field name in both outfile
out_kbar_name=str(sys.argv[11]) # kbar field name in both outfile
out_lat_name=bathy_inlat  # lat field name in GEBCO outfile
out_lon_name=bathy_inlon  # lon field name in GEBCO outfile

npy_hrms_pre=str(sys.argv[12]) # npy hrms Output file prename
npy_kbar_pre=str(sys.argv[13]) # npy kbar Output file prename

eas_bathy=str(sys.argv[14])       # Path/Name of the EAS system bathymetry file
eas_lat=str(sys.argv[15])         # eas lat field name in the mesh mask file
eas_lon=str(sys.argv[16])         # eas lon field name in the mesh mask file
eas_Bathymetry=str(sys.argv[17])  # eas Bathymetry field name in the bathy file

eas_mesh=str(sys.argv[18])  # Path/Name of the EAS system mesh mask file
eas_tmask=str(sys.argv[19])  # eas tmask field name in the mesh mask file

outfile_eas=workdir+'/'+str(sys.argv[20]) # EAS-grid Output file name

# --- SET PARAMETERS

# set FIGURE info
figdir  = workdir+'/plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

# set BOXES info for GEBCO 1/30'  
resol     = 120 # Npt/deg = 1'
boxdim    = 5   # deg

# set FLAG
flag_calc_all  = False
flag_gebco_vars_save = False
flag_regrid    = True
flag_eas_vars_save = True
flag_hrms_plot = True
flag_Kbar_plot = True


#### ------- DO NOT CHANGE below this line ------- #####

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# Read GEBCO bathymetry file
print ('Reading: [%s]' %fn)
ncf = Dataset(fn,mode='r')
bathy = ncf.variables[bathy_inname][:]       ; bathy = np.squeeze(bathy)
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
if flag_calc_all:
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

   # Med subdivision in 40 subdomains
   DOMY = [0,480,720,961,1200,NY-1] 
   DOMX = [0,840,1680,2520,3360,4200,5040,5880,NX-1]

   # +--+--+--+--+--+--+--+--+
   # |  |  |  |  |  |  |  |  |
   # +--+--+--+--+--+--+--+--+
   # |  |  |  |  |  |  |  |  |
   # +--+--+--+--+--+--+--+--+
   # |  |  |  |  |  |  |  |  |
   # +--+--+--+--+--+--+--+--+
   # |  |  |  |  |  |  |  |  |
   # |  |  |  |  |  |  |  |  |
   # +--+--+--+--+--+--+--+--+


   # Merge subdomains
   for bx in range(0,sub_bx_num+1):
       for by in (0,2,3,4): #range(0,sub_by_num+1):
           subdom = str(by)+str(bx)
           print('I am merging the subdomain: ',subdom)

           tmp_h = np.load(workdir+'/'+npy_hrms_pre+'dom'+subdom+'.npy')
           tmp_k = np.load(workdir+'/'+npy_kbar_pre+'dom'+subdom+'.npy')
           try:
              hrms_f[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]] = tmp_h[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]] 
              kbar_f[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]] = tmp_k[DOMY[by]:DOMY[by+1],DOMX[bx]:DOMX[bx+1]]
           except:
              hrms_f[DOMY[by]:DOMY[by+2],DOMX[bx]:DOMX[bx+1]] = tmp_h[DOMY[by]:DOMY[by+2],DOMX[bx]:DOMX[bx+1]]
              kbar_f[DOMY[by]:DOMY[by+2],DOMX[bx]:DOMX[bx+1]] = tmp_k[DOMY[by]:DOMY[by+2],DOMX[bx]:DOMX[bx+1]]

   # apply mask to the final field (to be sure...)
   hrms_f = np.where(msk==1, hrms_f, -9999.)
   kbar_f = np.where(msk==1, kbar_f, -9999.)


#   ---------------------
#   | write netCDF file |
#   --------------------

if flag_gebco_vars_save :
   # open netCDF file to write
   nco  = Dataset(outfile_gebco,mode='w',format='NETCDF4')

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
   print ("Saving: [%s]" %outfile_gebco)


#   ----------------------------------------------
#  | Regrid from GEBCO (1/120) to EAS (1/24) grid |
#   ----------------------------------------------
if flag_regrid:
   # Read new grid structure from EAS mesh_mask and bathy_meter files

   # Bathy file
   filebathy = eas_bathy
   print ('Reading: [%s]' %filebathy)
   ncf = Dataset(filebathy,mode='r')
   nav_lat = ncf.variables[eas_lat][:]   ; nav_lat = np.squeeze(nav_lat)        # Y-axis
   nav_lon = ncf.variables[eas_lon][:]   ; nav_lon = np.squeeze(nav_lon)        # X-axis
   bathy = ncf.variables[eas_Bathymetry][:]   ; bathy = np.squeeze(bathy)
   ncf.close()
  
   # Mes mask file
   filemesh = eas_mesh
   print ('Reading: [%s]' %filemesh)
   ncf = Dataset(filemesh,mode='r')
   tmask = ncf.variables[eas_tmask][0,0,:,:]   ; tmask = np.squeeze(tmask)
   ncf.close()

   # Defn of new (EAS) grid
   x=nav_lon[0] # Array 1D longitudes
   y=[ el[0] for el in nav_lat] #  Array 1D latitudes
   y=np.asarray(y)


   # Read the fields to be interpolated
   ncf = Dataset(outfile_gebco,mode='r')
   hrms  = ncf.variables[out_hrms_name][:]      ; hrms= np.squeeze(hrms)
   kbar  = ncf.variables[out_kbar_name][:]      ; kbar= np.squeeze(kbar)
   ncf.close()

   # Apply the SeaOverLand before interpolating

   for tobeSOL in (hrms,kbar):

       threshhk = -9999
       maskhk = globals()[tobeSOL] == threshhk
       fieldhk_ma = np.ma.masked_where(maskhk,globals()[tobeSOL])
       input_matrix=fieldhk_ma

       print ('I am applying SeOverLand to the following field: ',tobeSOL)

       nloop=napp_bathy
       infill_value = input_matrix.fill_value
       output_matrix = ma.copy(input_matrix)  # using ma.copy to prevent future field modifications

       if np.sum(output_matrix.mask) == 0:  # nothing to fill
          print('WARNING. Field does not have any land points or wrong field type. Exiting.', file=sys.stderr)
       else:
       # iterations loop
          for loop in range(nloop):
             # Create a nD x 8 matrix in which, last dimension fixed, the other dimensions
             #  contains values that are shifted in one of the 8 possible directions
             # of the last two axes compared to the original matrix
             shift_matrix = ma.array(np.empty(shape=((8,) + output_matrix.shape)),
                                     mask=True, fill_value=infill_value, dtype=float)
             # up shift
             shift_matrix[0, ..., : -1,:] = output_matrix[..., 1:,:]
             # down shift
             shift_matrix[1, ..., 1:, :] = output_matrix[..., : -1, :]
             # left shift
             shift_matrix[2, ..., :, : -1] = output_matrix[..., :, 1:]
             # right shift
             shift_matrix[3, ..., :, 1:] = output_matrix[..., :, : -1]
             # up-left shift
             shift_matrix[4, ..., : -1, : -1] = output_matrix[..., 1:, 1:]
             # up-right shift
             shift_matrix[5, ..., : -1, 1:] = output_matrix[..., 1:, : -1]
             # down-left shift
             shift_matrix[6, ..., 1:, : -1] = output_matrix[..., : -1, 1:]
             # down-right shift
             shift_matrix[7, ..., 1:, 1:] = output_matrix[..., : -1, : -1]
             # Mediate the shift matrix among the third dimension
             mean_matrix = ma.mean(shift_matrix, 0)
             # Replace input masked points with new ones belonging to the mean matrix
             output_matrix = ma.array(np.where(mean_matrix.mask + output_matrix.mask, mean_matrix, output_matrix),
                                      mask=mean_matrix.mask, fill_value=infill_value, dtype=float)
             output_matrix = ma.masked_where(mean_matrix.mask, output_matrix)
             if np.sum(output_matrix.mask) == 0:  # nothing more to flood
                 print('WARNING. Field does not have anymore land points,', str(loop + 1),
                       'steps were sufficient to flood it completely.', file=sys.stderr)
                 break
       # Sobstitute the SOL field
       globals()[tobeSOL]=output_matrix
       print('..Done!')

   # Set the interpolation
   print ('I am going to interpolate hrms and kbar fields from GEBCO to EAS grid..')
   h_field2interp=hrms
   k_field2interp=kbar

   h_field2interp=np.reshape(h_field2interp,len(lon)*len(lat))
   k_field2interp=np.reshape(k_field2interp,len(lon)*len(lat))

   # Interpolation from the old to the new grid 
   h_new_tmp = interpolate.interp2d(lon,lat,h_field2interp)
   k_new_tmp = interpolate.interp2d(lon,lat,k_field2interp)

   hrms_eas = h_new_tmp(x,y)
   kbar_eas = k_new_tmp(x,y)
   print('..Done!')

   # Mask the new field with the correct tmask
   print ('I am going to mask the fields with the correct tmask..')  
 
   # Fill nans with 0 
   hrms_eas=np.nan_to_num(hrms_eas)
   hrms_eas=np.where(hrms_eas==infill_value,0,hrms_eas)

   kbar_eas=np.nan_to_num(kbar_eas)
   kbar_eas=np.where(kbar_eas==infill_value,0,kbar_eas)

   # Apply the Mask 
   hrms_eas[time_output,:,:]=hrms_eas[time_output,:,:]*tmask[:,:]
   kbar_eas[time_output,:,:]=kbar_eas[time_output,:,:]*tmask[:,:]
   print('..Done!')


   # Add the outfields in the outfile_eas
   if flag_eas_vars_save :

   nc2open=workdir+'/'+outfile_eas
   print ('I am going to open and modify the following file: ',nc2open)
   temp_file = NC.Dataset(nc2open,'r+')

   # Create the new fields to be appended in the input/output netCDF
   print ('I am going to add the new regridded hrms and kbar fields to the outfile template ',outfile_eas)
   new_eas_hrms=temp_file.createVariable(out_hrms_name,np.float64,(nav_lat,nav_lon))
   new_eas_hrms.units = 'm'
   new_eas_kbar=temp_file.createVariable(out_kbar_name,np.float64,(nav_lat,nav_lon))
   new_eas_kbar.units = 'm-1'

   # Add the fields
   new_eas_hrms[:]=hrms_eas[:]
   new_eas_kbar[:]=kbar_eas[:]

   temp_file.close()
   print ('Done')

#   -----------------------
#   |    make the plot    |
#   -----------------------

if flag_hrms_plot :

   # open netCDF file to read
   #ncf = Dataset(outfile,mode='r')
   #hrms_eas  = ncf.variables[out_hrms_name][:]      ; hrms= np.squeeze(hrms)
   #kbar_eas  = ncf.variables[out_kbar_name][:]      ; kbar= np.squeeze(kbar) 

   # Read bathymetry file
   filebathy = eas_bathy
   print ('Reading: [%s]' %filebathy)
   ncf = Dataset(filebathy,mode='r')
   nav_lat = ncf.variables[eas_lat][:]   ; nav_lat = np.squeeze(nav_lat)        # Y-axis
   nav_lon = ncf.variables[eas_lon][:]   ; nav_lon = np.squeeze(nav_lon)        # X-axis
   bathy = ncf.variables[eas_Bathymetry][:]   ; bathy = np.squeeze(bathy)
   ncf.close()

   # Area
   # Mask invalid values
   #nav_lon[nav_lon==0]=np.nan
   # Mask variables
   nav_lat = np.ma.masked_invalid(nav_lat)
   lat_min=np.min(nav_lat[:,0])
   lat_max=np.max(nav_lat[:,0])
   lon_min=np.min(nav_lon[0,:])
   lon_max=np.max(nav_lon[0,:])
   nav_lon = np.ma.masked_invalid(nav_lon)   
   
   VAR = hrms_eas
   VARunit = '[m]'   
   VAR = np.ma.masked_invalid(VAR)

   figname = figdir +'map_hrms.png'
   figtitle = r'$h_{rms} = \sqrt{\frac{1}{4A\pi^2} \int\int |\hat{h}|^2 dk dl}$'
   cmap        = plt.cm.gist_heat_r # Colormap
   [cmin,cmax] = [0,600]            # color min and max values
      
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   #fig   = plt.contourf(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   pcf   = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0], colors='black',linewidth=0.3) # tmask[0,:,:]
   plt.title( figtitle, fontsize='18')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
   
   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')

if flag_Kbar_plot :
   VAR = kbar_eas
   VARunit = r'[$m^{-1}]'
   VAR = np.ma.masked_invalid(VAR)

   figname = figdir +'map_Kbar.png'
   figtitle = r'$\overline{K} = \frac{1}{4A\pi^2h_{rms}^2} \int\int (k^2+l^2) |\hat{h}|^2 dk dl$'
   cmap        = plt.cm.gist_heat_r   # Colormap
   [cmin,cmax] = [0.,1.e-4]       # color min and max values

   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, vmin=cmin, vmax=cmax)
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0], colors='black',linewidth=0.3) # tmask[0,:,:]
   plt.title( figtitle, fontsize='18')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
   cbar.formatter.set_powerlimits((0, 0))

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_inches='tight')
   plt.close('all')

   ncf.close()
#
