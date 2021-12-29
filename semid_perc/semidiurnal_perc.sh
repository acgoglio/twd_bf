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
WORKDIR=$1          # Work directory
PERC_INFILE=$2      # Input path/file with Amplitudea and Phases of all the tidal components
AMPPERC_OUTFILE=$3  # Name of the outfile name storing the ampperc field
AMPPERC_FIELD=$4    # Semidiurnal tidal amplitude percentage field

###############################
# Check the inputs and Move to the workdir
if [[ -d $WORKDIR ]]; then
   cd $WORKDIR
else
   echo "ERROR: WORKDIR=$1 NOT Found!"
   exit
fi

# Check the input file
if [[ ! -e $PERC_INFILE ]]; then
   echo "ERROR: Input file: $PERC_INFILE NOT FOUND! ..Why?"
   exit
fi

# Intermediate file handling
INTERMEDIATE_FILE="$WORKDIR/Amp_perc_tmp.nc"
if [[ -e $INTERMEDIATE_FILE ]]; then
   rm -v $INTERMEDIATE_FILE
fi


# Compute the Semidiurnal tidal amplitude percentage field
cdo expr,"ST_Amp=M2_Amp+S2_Amp+N2_Amp+K2_Amp;TOT_Amp=M2_Amp+S2_Amp+N2_Amp+K2_Amp+K1_Amp+O1_Amp+P1_Amp+Q1_Amp;$AMPPERC_FIELD=100.0*ST_Amp/TOT_Amp" $PERC_INFILE $INTERMEDIATE_FILE

# Set the field units
# ncatted -a units,'ST_perc',c,c,'%' $INTERMEDIATE_FILE

# Set missing values to zeros
cdo setmisstoc,0 $INTERMEDIATE_FILE $AMPPERC_OUTFILE

#echo "Do some cleaning of intermediate files.."
rm $INTERMEDIATE_FILE

echo "All done!"
