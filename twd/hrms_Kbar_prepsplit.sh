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
MESH_FILE=$2   # EAS mesh mask path and file name
MESH_TMASK=$3  # tmask field name
HK_EAS_OUTFILE=$4  # Name of the template file where to store hrms and kbar on the EAS grid
TASK_FLAG=$5       # Flag to select between tasks: B=build template intermediate file ; S=split intermediate file
HRMS_OUTFILE=$6    # Outfile storing hrms field
KBAR_OUTFILE=$7    # Outfile storing kbar field
HRMS_OUTFIELD=$8   # hrms field
KBAR_OUTFIELD=$9   # kbar field

###############################
# Check the inputs and Move to the workdir
if [[ -d $WORKDIR ]]; then
   cd $WORKDIR
else
   echo "ERROR: WORKDIR=$1 NOT Found!"
   exit
fi

if [[ ! -e ${MESH_FILE} ]]; then
   echo "ERROR: MESH_FILE=$2 NOT Found!"
   exit
fi
INTERMEDIATE_FILE=$WORKDIR/$HK_EAS_OUTFILE
if [[ -e $INTERMEDIATE_FILE ]]; then
   rm -v $INTERMEDIATE_FILE
fi
if [[ $TASK_FLAG == 'B' ]]; then 
   echo "I am going to build the 2ND INTERMEDIATE OUTPUT (HRMS and KBAR FIELDS on EAS GRID).."
elif [[ $TASK_FLAG == 'S' ]]; then
   echo "I am going to split the intermediate file into the 2 outfiles ).."
else
   echo "ERROR: Wrong FLAG: use B to build template intermediate file or S to split intermediate file "
   exit
fi

# CASE: build the 2ND INTERMEDIATE OUTPUT (HRMS and KBAR FIELDS on EAS GRID)
if [[ $TASK_FLAG == 'B' ]]; then
   # build the intermediate file storing hrms and kbar on the eas grid from the mesh mask
   echo "I am building the file template.."
   cdo select,name=$MESH_TMASK,timestep=1 $MESH_FILE $INTERMEDIATE_FILE 
   echo "Done!"
fi

# CASE: split the intermediate file into the 2 outfiles
if [[ $TASK_FLAG == 'S' ]]; then
   # split the intermediate file into the 2 outfiles
   echo "I am splitting the intermediate file.."
   cdo select,name=$HRMS_OUTFIELD $INTERMEDIATE_FILE $HRMS_OUTFILE
   cdo select,name=$KBAR_OUTFIELD $INTERMEDIATE_FILE $KBAR_OUTFILE
   echo "Done!"
fi

echo "Do some cleaning of intermediate files.."
#rm $INTERMEDIATE_FILE

echo "All done!"
