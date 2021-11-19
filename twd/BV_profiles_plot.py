"""
Written by AC Goglio , Nov 2021
"""

print ('\n $$$--- Plot frofiles of Brunt-Vaisala Frequency ---$$$\n')

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

for season in ('yearly','DJF','MAM','JJA','SON'):
   for where in ('3035','3540','4045'):
      name_file=season+'_mean_bn2.nc_'+where+'.nc_S.nc'
      # --- READ ARGS
      
      workdir='/work/oda/ag15419/tmp/BV_N2_profiles/' #str(sys.argv[1])   # Work directory
      eas_mesh='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/mesh_mask.nc' #str(sys.argv[2])  # Path/Name of the EAS system mesh mask file
      nav_lev='nav_lev' #str(sys.argv[3])   # eas nav_lev field name in the mesh mask file
      
      outfile=workdir+'/'+name_file #str(sys.argv[2])   # Input file name
      outfield='bn2' #str(sys.argv[7])   # Name of the n2 field in the input file
      
      #where=str(sys.argv[3])  # Name of the box
      #season=str(sys.argv[4]) # Period: year or season
      
      
      ##############
      
      # --- SET PARAMETERS
      flag_plot_shap=1
      
      # set FIGURE info
      figdir  = workdir+'/plots/'
      if not(os.path.isdir(figdir)) :
         print('Creating: [%s]' %figdir)
         os.makedirs(figdir)
      
      if flag_plot_shap :
         # Read file
         fn = outfile
         print ('Reading: [%s]' %fn)
         ncf = Dataset(fn,mode='r')
         bn2 = ncf.variables[outfield][:]   
         bn2 =np.squeeze(bn2)
         bn=np.sqrt(bn2)
         ncf.close()
      
         mesh_nc = Dataset(eas_mesh,mode='r')
         print ('Reading: [%s]' %mesh_nc)
         lev_m = mesh_nc.variables[nav_lev][:]
         lev_m=-1.0*lev_m
         mesh_nc.close()
      
      
         # --- PLOT
         VAR=bn
         VARunit = r'[$s^{-1}$]'
         #VAR = np.ma.masked_invalid(VAR)
         figname = figdir +'profile_'+season+'_'+where+'.png'
         figtitle = 'N Profile '+season+'_'+where
      
         print('... make the plot ...')
      
         plt.figure(figsize=(7,10))
         plt.rc('font', size=16)
         plt.title (figtitle)
         plt.plot(VAR,lev_m,color='#d62728',linewidth=3,label = 'N Profile')
         plt.grid ()
         plt.ylabel ('Depth [m]')
         plt.xlabel ('N '+VARunit)
         plt.legend()
         plt.xlim(0e-3,20e-3)
         plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
         plt.ylim(-2000,50)
      
         print ('Saving: [%s]' % figname)
         plt.savefig(figname)
         plt.clf()
      
      #
      ## --- END
      #
