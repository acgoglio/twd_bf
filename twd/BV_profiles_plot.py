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
season_Fullname=['','Winter','Spring','Summer','Autumn']
where_num=['366','382','398']

seasons=['yearly','DJF','MAM','JJA','SON']
wheres=['3035','3540','4045']

for idx_season,season in enumerate(seasons):
   for idx_where,where in enumerate(wheres):
      
      name_file=season+'_mean_bn2.nc_'+where+'.nc_S.nc'

      # --- READ ARGS
      
      workdir='/work/oda/ag15419/tmp/BV_N2_profiles/' #str(sys.argv[1])   # Work directory
      eas_mesh='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/mesh_mask.nc' #str(sys.argv[2])  # Path/Name of the EAS system mesh mask file
      nav_lev='nav_lev' #str(sys.argv[3])   # eas nav_lev field name in the mesh mask file
      
      outfile=workdir+'/'+name_file #str(sys.argv[2])   # Input file name
      outfield='bn2' #str(sys.argv[7])   # Name of the n2 field in the input file
      
      #where=str(sys.argv[3])  # Name of the box
      #season=str(sys.argv[4]) # Period: year or season
      
      if season != 'yearly':
         obspath=workdir+'/OBS_prof/'
         obsname='SDC_ATL_N2_'+where_num[idx_where]+'_'+season_Fullname[idx_season]+'.nc'
         obs_field='BVF'
         obs_depth='depth'
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
      
         if season != 'yearly':
            obs_nc = Dataset(obspath+'/'+obsname,mode='r')
            print ('Reading: [%s]' %obs_nc)
            bn2_obs = obs_nc.variables[obs_field][:]
            bn2_obs = np.squeeze(bn2_obs)
            bn_obs = np.sqrt(bn2_obs)
            depth_obs = obs_nc.variables[obs_depth][:]   
            depth_obs = -1.0*np.squeeze(depth_obs)
            obs_nc.close()

         # --- PLOT
         VAR=bn
         VARunit = r'[$s^{-1}$]'
         #VAR = np.ma.masked_invalid(VAR)
         figname = figdir +'profile_'+season+'_'+where+'.png'
         if season != 'yearly':
            figtitle = 'N Profile '+season_Fullname[idx_season]+' '+where
         else:
            figtitle = 'N Profile '+'yearly mean'+' '+where      

         print('... make the plot ...')
      
         plt.figure(figsize=(7,10))
         plt.rc('font', size=16)
         plt.title (figtitle)

         if season != 'yearly':
            plt.plot(VAR,lev_m,color='#d62728',linewidth=3,label = 'Model')
            plt.plot(bn_obs,depth_obs,color='#1f77b4',linewidth=3,label = 'Obs') 
         else:
            plt.plot(VAR,lev_m,color='#d62728',marker='o',label = 'N Profile')

         plt.grid ()
         plt.ylabel ('Depth [m]')
         plt.xlabel ('N '+VARunit)
         plt.legend()

         if season != 'yearly':
            plt.xlim(0e-3,20e-3)

         plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
         
         if season != 'yearly':
            plt.ylim(-2000,50)

         print ('Saving: [%s]' % figname)
         plt.savefig(figname)
         plt.clf()
      
      #
      ## --- END
      #
