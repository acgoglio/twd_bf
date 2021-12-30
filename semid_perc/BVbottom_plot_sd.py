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
eas_gphif='gphif'

order=int(sys.argv[13]) # Order of the Shapiro filter
napp=int(sys.argv[14]) # Number of Shapiro filter applications
scheme=int(sys.argv[15]) # type of boundary scheme to use ( only option 1 implemented = No change at wall, constant order )  

hk_outfile=workdir+'/'+str(sys.argv[16]) # File storing hrms and kbar for the TWD coeff plot
out_hrms_name=str(sys.argv[17]) # hrms field name in both outfile
out_kbar_name=str(sys.argv[18]) # kbar field name in both outfile

# Name of the semidiurnal tidal amplitude percentage input file 
semid_file=str(sys.argv[19])
# Name of the semidiurnal tidal amplitude percentage input field
semid_field=str(sys.argv[20])

# --- SET PARAMETERS

# set FIGURE info
figdir  = workdir+'/plots/'
if not(os.path.isdir(figdir)) :
   print('Creating: [%s]' %figdir)
   os.makedirs(figdir)

# set FLAG
flag_calc_bnbot = False
flag_save_bnbot = False
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
gphif   = mesh.variables[eas_gphif][:]     ; gphif   = np.squeeze(gphif)
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
   #bn=ma.masked_array(bn,mask = tmask,fill_value=np.nan)
   bnb = np.empty([NY,NX])
   for i in range(NX) :
      for j in range(NY) :
         bmk = mbathy[j,i]-1
         #print ('bmk',bmk)
         if bmk>=0 :
            bnb[j,i] = bn[bmk,j,i]
         else: 
            #print ('mbathy[j,i]',mbathy[j,i])
            bnb[j,i] = -1 #np.nan


# --------- Apply SeaOver Land before smoothing
   # Sea over land
   threshbnb = -1
   maskbnb = bnb == threshbnb
   fieldbnb_ma = np.ma.masked_where(maskbnb,bnb)
   input_matrix=fieldbnb_ma
   nloop=napp
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
   bnb=output_matrix


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
if napp != 0 and flag_calc_bnbot :
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

if flag_plot_shap :
   # Read file
   fn = outfile
   print ('Reading: [%s]' %fn)
   ncf = Dataset(fn,mode='r')
   bnbot = ncf.variables[outfield][:]   
   print ('Done')
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
   
   # --- PLOT Bottom Brunt Vaisal freq.
   VAR = bnbot
   VARunit = r'[$s^{-1}$]'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_bottomBV.png'
   figtitle = r'$N_{bottom}$'
   cmap        = cm.get_cmap('viridis')
   [cmin,cmax] = [1.e-5,1.e-1]
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3) 
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.pcolor(x,y,VAR, cmap=cmap, norm=colors.LogNorm(vmin=cmin, vmax=cmax) ) # [0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1.0]
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3) 
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
   
   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')


   ################################

   # OPEN FILES AND READ NEEDED FIELDS

   # Open hk file path/name
   nc2open=hk_outfile
   print ('I am going to open and plot the following file: ',nc2open)
   bathy_infield = Dataset(nc2open,'r')
   hrms=bathy_infield.variables[out_hrms_name][:]
   kbar=bathy_infield.variables[out_kbar_name][:]
   bathy_infield.close()

   # Read semidiurnal percentages
   nc2open=workdir+'/'+semid_file
   print ('Input file = ',nc2open)
   semid_mod = Dataset(nc2open,'r')
   ST_perc=semid_mod.variables[semid_field][:]

   # DEFINE NEEDED PARAMATER
   Omega = 1.405189*0.0001 # s-1 M2 angular frequency 
   Omega_k1 = 7.292117*0.00001 # s-1 K1 angular frequency

   f_cor = 2*7.2921*0.00001*np.sin(gphif*(2*np.pi)/360.)  # s-1 Coriolis param

   # COMPIUTE THE VALUES TO BE PLOTTED
   W_fun = np.where((bnbot*bnbot-Omega*Omega)*(Omega*Omega-f_cor*f_cor)>0,np.sqrt((bnbot*bnbot-Omega*Omega)*(Omega*Omega-f_cor*f_cor))/(Omega*bnbot),0)
   W_fun_k1 = np.where((bnbot*bnbot-Omega_k1*Omega_k1)*(Omega_k1*Omega_k1-f_cor*f_cor)>0,np.sqrt((bnbot*bnbot-Omega_k1*Omega_k1)*(Omega_k1*Omega_k1-f_cor*f_cor))/(Omega_k1*bnbot),0)

   sign_W_fun = np.sign((bnbot*bnbot-Omega*Omega)*(Omega*Omega-f_cor*f_cor))
   sign1_W_fun = np.sign(bnbot*bnbot-Omega*Omega)
   sign2_W_fun = np.sign(Omega*Omega-f_cor*f_cor)
   sign_W_fun_k1 = np.sign((bnbot*bnbot-Omega_k1*Omega_k1)*(Omega_k1*Omega_k1-f_cor*f_cor))

   TWD_coeff = 0.5*(hrms*hrms)*kbar*W_fun*bnbot
   TWD_coeff_k1 = 0.5*(hrms*hrms)*kbar*W_fun_k1*bnbot

   SDandD_W_fun=(W_fun*(ST_perc/100.0))+(W_fun_k1*((100.0-ST_perc)/100.0))

   sign_W_fun_SDandD = np.sign(SDandD_W_fun) 

   TWD_coeff_SDandD = 0.5*(hrms*hrms)*kbar*SDandD_W_fun*bnbot

   # DEFINE THE LEVELS
   W_fun_levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

   [cmin,cmax] = [0.0000001,0.01]
   [cmin_diff,cmax_diff] = [0.00000001,0.001]

   TWD_levels=[0.0000001,0.000001,0.00001,0.0001,0.001,0.01]
   TWD_levels_diff=[0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001]

   TWD_diff_levels=[-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9]
   #[-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #[-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5]

   ################################
   # --- PLOT M2 Weighting Function

   VAR = W_fun
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_weightingF.png'
   figtitle = 'Weighting Function (M2 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=W_fun_levels,cmap=cmap,extend='max') 
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOT K1 Weighting Function

   VAR = W_fun_k1
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_weightingF_k1.png'
   figtitle = 'Weighting Function (K1 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=W_fun_levels,cmap=cmap,extend='max')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')


   # --- PLOT Sign of the Weighting Function

   VAR = sign_W_fun
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_sign_weightingF.png'
   figtitle = 'Sign of the Weighting Function (M2 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=[-1.0,0.0,1.0],cmap=cmap )
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOT Sign of the K1 Weighting Function

   VAR = sign_W_fun_k1
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_sign_weightingF_k1.png'
   figtitle = 'Sign of the Weighting Function (K1 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=[-1.0,0.0,1.0],cmap=cmap )
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOT Sign 1 of the Weighting Function

   VAR = sign1_W_fun
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_sign_weightingF_bnbot-Omega.png'
   figtitle = 'Sign of bnbot*bnbot-Omega*Omega (M2 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   fig = m.contourf(x,y,VAR,levels=[-1.0,0.0,1.0],cmap=cmap )
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOT Sign 2 of the Weighting Function

   VAR = sign2_W_fun
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_sign_weightingF_Omega-f.png'
   figtitle = 'Sign of Omega*Omega-f_cor*f_cor (M2 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=[-1.0,0.0,1.0],cmap=cmap )
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   ##############################################

   # PLOT TWD semidiurnal coeff
   #
   VAR = TWD_coeff
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_TWDcoeff.png'
   figtitle = 'TWD coeff (M2 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_levels,cmap=cmap,norm=colors.LogNorm(vmin=cmin,vmax=cmax),extend='max')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOT TWD Coeff only for H>500m
   #
   VAR = TWD_coeff
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_TWDcoeff_500.png'
   figtitle = 'TWD coeff (M2 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_levels,cmap=cmap,norm=colors.LogNorm(vmin=cmin,vmax=cmax),extend='max')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOT TWD Coeff for K1 tidal component only for H>500m
   #
   VAR = TWD_coeff_k1
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_TWDcoeff_500_k1.png'
   figtitle = 'TWD coeff (K1 tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_levels,cmap=cmap,norm=colors.LogNorm(vmin=cmin,vmax=cmax),extend='max')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOTDiff  TWD Coeff for diurnal-semidiurnal tidal component only for H>500m
   #
   VAR = TWD_coeff_k1-TWD_coeff
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'diff_TWDcoeff_500_M2-K1.png'
   figtitle = 'Diff: TWD coeff diurnal - TWD coeff semidiurnal'
   cmap        = cm.get_cmap('bwr')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_diff_levels,cmap=cmap,extend='both')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOTDiff  TWD Coeff for diurnal-semidiurnal tidal component only for H>500m LOG scale
   #
   VAR = TWD_coeff_k1-TWD_coeff
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'difflog_TWDcoeff_500_M2-K1.png'
   figtitle = 'Diff: TWD coeff diurnal - TWD coeff semidiurnal'
   cmap        = cm.get_cmap('Reds')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_levels_diff,cmap=cmap,norm=colors.LogNorm(vmin=cmin_diff,vmax=cmax_diff),extend='max')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   ######### COMBINED SEMIDIURNAL + DIURNAL CONTRIBUTIONS ##############

   # ---- PLOT the new Function for diurnal+semidiurnal tidal components
   # Wheighted with amplitude percentages
   # only for H>500m

   # Define and mask the field to be plotted
   VAR = SDandD_W_fun
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_weightingF_m2k1.png'
   figtitle = 'Weighting Function (M2+K1 tidal components)'
   cmap        = cm.get_cmap('bone')
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=W_fun_levels,cmap=cmap,extend='max' )
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # ---- PLOT the DIFF between new Function for diurnal+semidiurnal tidal components
   # Wheighted with amplitude percentages and old Function for semidiurnal tidal component
   # only for H>500m

   # Define and mask the field to be plotted
   VAR =SDandD_W_fun-W_fun
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'diffmap_weightingF_m2k1-m2.png'
   figtitle = 'Weighting Function M2+K1 - Weighting Function M2'
   cmap        = cm.get_cmap('bwr')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=TWD_diff_levels,cmap=cmap,extend='both' )
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')


   # --- PLOT Sign of the new Weighting Function for diurnal+semidiurnal tidal components
   # Wheighted with amplitude percentages

   VAR = sign_W_fun_SDandD
   VARunit = ''
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   k = 0
   figname = figdir +'map_sign_weightingF_m2k1.png'
   figtitle = 'Sign of the Weighting Function (M2_K1 tidal components)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig = m.contourf(x,y,VAR,levels=[-1.0,0.0,1.0],cmap=cmap )
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,15.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')


   # ---- PLOT the new TWD coeff for diurnal+semidiurnal tidal components
   # Wheighted with amplitude percentages
   # only for H>500m

   VAR = TWD_coeff_SDandD
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'map_TWDcoeff_500_k1m2.png'
   figtitle = 'TWD coeff (semidiurnal+diurnal tidal component)'
   cmap        = cm.get_cmap('jet')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_levels,cmap=cmap,norm=colors.LogNorm(vmin=cmin,vmax=cmax),extend='both')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')
 
   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOTDiff  NEW TWD Coeff for semidiurnal+diurnal tidal component -
   # - OLD Coeff for semidiurnal tidal componentonly 
   # only for H>500m

   VAR = TWD_coeff_SDandD-TWD_coeff
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'diffmap_TWDcoeff_500_k1m2.png'
   figtitle = 'Diff: TWD coeff diurnal+semidiurnal - TWD coeff semidiurnal'
   cmap        = cm.get_cmap('bwr')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_diff_levels,cmap=cmap,extend='both')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='both')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

   # --- PLOTDiff  TWD Coeff for diurnal+semidiurnal - semidiurnal tidal component only for H>500m LOG scale
   #
   VAR = TWD_coeff-TWD_coeff_SDandD
   VARunit = 'm/s'
   VAR = np.ma.masked_invalid(VAR)
   VAR=VAR*tmask[0,:,:]

   figname = figdir +'difflog_TWDcoeff_M2-M2K1.png'
   figtitle = 'Diff: TWD coeff semidiurnal - TWD coeff diurnal+semidiurnal'
   cmap        = cm.get_cmap('Reds')
   print('... make the plot ...')
   plt.figure()
   plt.rcParams['lines.linewidth'] = 0.3
   m = Basemap(projection='mill',llcrnrlat=lat_min,urcrnrlat=lat_max,llcrnrlon=lon_min,urcrnrlon=lon_max,resolution='i')
   m.drawparallels(np.arange(30., 46., 5), labels=[1,0,0,0], fontsize=6,linewidth=0.3)
   m.drawmeridians(np.arange(-20., 40., 10), labels=[0,0,0,1], fontsize=6,linewidth=0.3)
   x, y = m(nav_lon, nav_lat)
   fig=plt.contourf(x,y,VAR,levels=TWD_levels_diff,cmap=cmap,norm=colors.LogNorm(vmin=cmin_diff,vmax=cmax_diff),extend='max')
   plt.contourf(x,y,VAR,levels=[-1000,0.0],colors='white')
   pcf  = plt.contourf(x,y,bathy, levels=[0.000,500.0], colors='dimgray')
   pc    = plt.contour(x,y,bathy, levels=[15.0,500], colors='black',linewidth=0.3)
   plt.title( figtitle, fontsize='16')
   cbar = m.colorbar(fig,'bottom', size='10%', pad='10%', extend='max')
   cbar.set_label(VARunit,fontsize='14')
   cbar.ax.tick_params(labelsize='12')

   print ('Saving: [%s]' % figname)
   plt.savefig(figname, dpi=500, bbox_Nptshes='tight')
   plt.close('all')

#
## --- END
#
