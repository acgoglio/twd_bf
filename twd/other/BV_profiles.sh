#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 29/10/2021
#
#set -u
set -e
#set -x 
###############################
# READ INPUTS
WORKDIR="/work/oda/ag15419/tmp/BV_N2_profiles"     # Work directory
EXP_PATH="/work/oda/ag15419/arc_link/simu_bn2/output/"    # Path of the experiment, e.g. /work/oda/ag15419/exp/eas6_v7/simu_bn2/output
EXP_NAME="simu_bn2"    # Name of the experiment, e.g. simu_bn2
FIELD_NAME="bn2"  # Name of the bn2 field, e.g. bn2
OUTFIELD="mean_bn2"    # Name of the rmsq field, e.g. bnbot
OUTFILE="bn2_prof.nc"     # Name of the outfile storing the rmsq field, e.g. BNBOT.nc
###############################
# Check the inputs and Move to the workdir
if [[ -d $WORKDIR ]]; then
   cd $WORKDIR
else
   echo "ERROR: WORKDIR=$WORKDIR NOT Found!"
   exit
fi

if [[ ! -d ${EXP_PATH} ]]; then
   echo "ERROR: EXP_PATH=$EXP_PATH NOT Found!"
   exit
else
   ALL_EXP="${EXP_PATH}/*/${EXP_NAME}_1d_*_grid_W.nc"
   for TOCAT in $( ls ${ALL_EXP} ); do 
      echo "Found NEMO outfile: $TOCAT " 
   done
fi

# Select the field
#IDX_TOCAT=0
#echo "Select the field ${FIELD_NAME}"
#ALL_EXP="${EXP_PATH}/*/${EXP_NAME}_1d_*_grid_W.nc"
#for TOCAT in $( ls ${ALL_EXP} ); do 
#    IDX_TOCAT=$(( $IDX_TOCAT + 1 ))
#    cdo select,name=${FIELD_NAME} $TOCAT selbn2_${IDX_TOCAT}.nc
#done 

# Cat all the files
#echo "Cat all the files"
#cdo mergetime selbn2_*.nc catted.nc
#rm selbn2_*.nc

# Compute the TIME mean 
#echo "Compute the time means:"
#echo "yearly mean.."
#cdo timmean catted.nc yearly_mean_bn2.nc
#echo "seasonal mean.."
#SEASON_MONTHS=("DJF" "MAM" "JJA" "SON")
#for M_IDX in ${SEASON_MONTHS[@]}; do
#           echo "Working on season: $M_IDX .."
#           MAP3D_SEASON_OUTFILE="${M_IDX}_mean_bn2.nc"
#           cdo selseason,${M_IDX} catted.nc ${M_IDX}_tmp.nc
#           cdo timmean ${M_IDX}_tmp.nc $MAP3D_SEASON_OUTFILE
#           #rm ${M_IDX}_tmp.nc
#done
#rm catted.nc

# Select the area and compute the SPACE mean 
#for TOCUT in $( ls *_mean_bn2.nc ); do
#    cdo sellonlatbox,-15,-10,30,35 $TOCUT ${TOCUT}_3035.nc
#    cdo sellonlatbox,-15,-10,35,40 $TOCUT ${TOCUT}_3540.nc
#    cdo sellonlatbox,-15,-10,40,45 $TOCUT ${TOCUT}_4045.nc
#done

# Compute the SPACE mean
echo "Compute the means over the boxes:"
for TOMEAN in $( ls *_mean_bn2.nc_*.nc ); do
   cdo fldmean $TOMEAN ${TOMEAN}_S.nc
   # rm $TOMEAN
done 

echo "All done!"
