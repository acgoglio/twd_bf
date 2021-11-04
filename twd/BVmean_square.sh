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
WORKDIR=$1     # Work directory
EXP_PATH=$2    # Path of the experiment, e.g. /work/oda/ag15419/exp/eas6_v7/simu_bn2/output
EXP_NAME=$3    # Name of the experiment, e.g. simu_bn2
FIELD_NAME=$4  # Name of the bn2 field, e.g. bn2
OUTFIELD=$5    # Name of the rmsq field, e.g. bnbot
OUTFILE=$6     # Name of the outfile storing the rmsq field, e.g. BNBOT.nc
###############################
# Check the inputs and Move to the workdir
if [[ -d $WORKDIR ]]; then
   cd $WORKDIR
else
   echo "ERROR: WORKDIR=$1 NOT Found!"
   exit
fi

if [[ ! -d ${EXP_PATH} ]]; then
   echo "ERROR: EXP_PATH=$2 NOT Found!"
   exit
else
   ALL_EXP="${EXP_PATH}/*/${EXP_NAME}_1d_*_grid_W.nc"
   for TOCAT in $( ls ${ALL_EXP} ); do 
      echo "Found NEMO outfile: $TOCAT "
   done
fi

# Cat all the files
echo "Cat all the files"
cdo mergetime ${ALL_EXP} catted.nc

## Select the field
#echo "Select the field"
#cdo select,name=${FIELD_NAME} catted.nc all_bn2.nc

# Compute the mean
echo "Compute the time mean"
cdo timmean catted.nc mean_bn2.nc

# Comute the rmsq
echo "Compute the sqrt"
EXPR_STRING="\'${OUTFIELD}=sqrt(${FIELD_NAME})\'"
cdo expr,$EXPR_STRING mean_bn2.nc ${OUTFILE}

echo "Do some cleaning of intermediate files.."
rm catted.nc 
rm mean_bn2.nc

echo "All done!"
