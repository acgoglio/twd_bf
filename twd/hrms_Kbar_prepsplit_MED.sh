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
HK_EAS_OUTFILE=$2  # Name of the template file where to store hrms and kbar on the EAS grid
HRMS_OUTFILE=$3    # Outfile storing hrms field
KBAR_OUTFILE=$4    # Outfile storing kbar field
HRMS_OUTFIELD=$5   # hrms field
KBAR_OUTFIELD=$6   # kbar field

###############################
# Check the inputs and Move to the workdir
if [[ -d $WORKDIR ]]; then
   cd $WORKDIR
else
   echo "ERROR: WORKDIR=$1 NOT Found!"
   exit
fi

INTERMEDIATE_FILE=$WORKDIR/$HK_EAS_OUTFILE
if [[ ! -e $INTERMEDIATE_FILE ]]; then
   echo "ERROR: The input intermedaite file $INTERMEDIATE_FILE DOES NOT EXIST..Why?!"
   exit
fi

# Split the intermediate file into the 2 outfiles
echo "I am splitting the intermediate file.."
cdo select,name=$HRMS_OUTFIELD $INTERMEDIATE_FILE $HRMS_OUTFILE
cdo select,name=$KBAR_OUTFIELD $INTERMEDIATE_FILE $KBAR_OUTFILE
echo "Done!"

#echo "Do some cleaning of intermediate files.."
#rm $INTERMEDIATE_FILE

echo "All done!"
