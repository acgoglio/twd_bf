#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
# Written: 17/12/2021
# Last modified: 17/12/2021

#set -u
set -e
#set -x 
########################
# Read ini file
source run_twd.ini

# Define the vars
SOURCE_DIR=$pwd
EXEC_DIR=$workdir

# Load the Env
module load anaconda 
source activate $conda_env

#######################
echo "----------------------------"
echo "Pre-processing $datetocompute" 

# Define the work, in and out directories
echo "SOURCE_DIR=$SOURCE_DIR"
echo "EXEC_DIR=$EXEC_DIR"
#
if [[ ! -d $EXEC_DIR ]]; then
   mkdir $EXEC_DIR
fi
#
for NEWDIR in in out work; do
   if [[ ! -d $EXEC_DIR/$NEWDIR ]]; then
      mkdir $EXEC_DIR/$NEWDIR
   fi
done
WORK_DIR="$EXEC_DIR/work"
OUT_DIR="$EXEC_DIR/out"

# Load the env (or define the needed vars)

# Define and move to the test-case source dir
NAME_SUBDIR="twd"
REP_SUBDIR=${SOURCE_DIR}/${NAME_SUBDIR}/
echo "The src are in: ${REP_SUBDIR}"

NAME_SUBDIR2="semid_perc"
REP_SUBDIR2=${SOURCE_DIR}/${NAME_SUBDIR2}/
echo "Other src are in: ${REP_SUBDIR2}"

# Activate (=1 or 2 or 3 for single steps or all) or deactivate (=0) the production of the outfiles
ROUGHNESS_FLAG=2
BRUNTV_FLAG=0
DIURNAL_FLAG=0

############################ ROUGHNESS HRMS and KBAR ###############################################
if [[ $ROUGHNESS_FLAG == 1 ]] || [[ $ROUGHNESS_FLAG == 'all' ]]; then
   # Run the code
   echo "----------------------------"
   echo "Running the roughness procedure..."
   
   # ----------------------------------
   # 1) Compute the roughness field as the difference between the GEBCO on MEd grid and EAS bathymetry
   #
   SCRIPT_TO_RUN="roughness_MEDgrid.py"
   #
   # Cp the oiginal file to the workdir 
   echo "INPUT=$gebco_bathy/$bathy_infile"
   echo "OUTPUT=$WORK_DIR/$bathy_outfile"
   
   cp $gebco_bathy/$bathy_infile $WORK_DIR/$bathy_outfile
   
   # Set the roughness_shapiro job parameters
   JOUT=stdout_roughness_%J
   JERR=stderr_roughness_%J
   JNAME="ROUGHNESS"
   JQUEUE="s_short"
   
   # Submit the roughness_shapiro job
   echo "Submitting the job: python ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $bathy_outfile $bathy_inname $bathy_outname $bathy_rough $bathy_inlat $bathy_inlon $eas_bathy_meter $eas_bathymetry_name $eas_lat_name $eas_lon_name $np2be0"
   
   bsub -K -J $JNAME -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE \
   "python ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $bathy_outfile $bathy_inname $bathy_outname $bathy_rough $bathy_inlat $bathy_inlon $eas_bathy_meter $eas_bathymetry_name $eas_lat_name $eas_lon_name $np2be0"
   
   echo "Done!"
   
   # Cp the outputs to the archive
   cp $WORK_DIR/$bathy_outfile $OUT_DIR/ 
fi   
if [[ $ROUGHNESS_FLAG == 2 ]] || [[ $ROUGHNESS_FLAG == 'all' ]]; then
   # ----------------------------------
   # 2) Compute h_rms and K_bar on the MED grid and plot the fields 
   #
   SCRIPT_TO_RUN="hrms_Kbar_MED.py"
   #
   echo "INPUT=$WORK_DIR/$bathy_outfile"
   echo "INTERM OUTPUT=$WORK_DIR/$interm_outfile"
   echo "OUTPUTS=$WORK_DIR/$hrms_outfile and $WORK_DIR/$kbar_outfile"
   
   # Set the roughness_shapiro job parameters
   JOUT=stdout_hk_%J
   JERR=stderr_hk_%J
   JNAME="HK_ROUGH"
   JQUEUE="s_medium"

   # Submit the hrms kbar jobs (N jobs for N Med subdomains..)
   echo "Working on the Whole Domain.."

   # Build the file template
   cp $WORK_DIR/$bathy_outfile $WORK_DIR/$interm_outfile

   # Run the job
   echo "Submitting the job: python ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $bathy_outfile $bathy_hrms $bathy_kbar $bathy_outname $bathy_rough $eas_lat_name $eas_lon_name $out_hrms_name $out_kbar_name $interm_outfile $resol $boxdim $order_hk $napp_hk $scheme_hk"

   bsub -K -J ${JNAME} -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE -M 80G \
   "nohup python ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $bathy_outfile $bathy_hrms $bathy_kbar $bathy_outname $bathy_rough $eas_lat_name $eas_lon_name $out_hrms_name $out_kbar_name $interm_outfile $resol $boxdim $order_hk $napp_hk $scheme_hk"

   echo "Done!"
   
   # Cp the outputs to the archive
   cp $WORK_DIR/$interm_outfile $OUT_DIR/
fi

if [[ $ROUGHNESS_FLAG == 3 ]] || [[ $ROUGHNESS_FLAG == 'all' ]]; then
   # ----------------------------------
   # 3) Split the outfile in the single output files storing respectively hrms and kbar fields
   #
   SCRIPT_TO_RUN="hrms_Kbar_prepsplit_MED.sh"
   #
   # Set the splitting job parameters
   JOUT=stdout_split_%J
   JERR=stderr_split_%J
   JNAME="SPLIT_HK"
   JQUEUE="s_short"

   # Submit the pre proc job
   echo "Submitting the job: sh ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $interm_outfile $hrms_outfile $kbar_outfile $out_hrms_name $out_kbar_name"

   bsub -K -J $JNAME -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE \
   "sh ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $interm_outfile $hrms_outfile $kbar_outfile $out_hrms_name $out_kbar_name"

   # ---
   # Store the outputs
   # Cp the outputs to the archive
   cp $WORK_DIR/$hrms_outfile $OUT_DIR/
   cp $WORK_DIR/$kbar_outfile $OUT_DIR/

fi

################################### BOTTOM BRUNT VAISALA ##########################################
# Name of the intermediat outfile
INT_OUTFILE="${bvf_name}.nc"

if [[ $BRUNTV_FLAG == 1 ]] || [[ $BRUNTV_FLAG == 'all' ]]; then
 
   # Run the code
   echo "----------------------------"
   echo "Running the bottom Brunt Vaisala procedure..."
   # ----------------------------------
   # 1) Compute the mean and rmsq of the squared brunt vaisala field in the NEMO outputs
   echo "INPUT (exp path and name) : $nontidal_exp_path $nontidal_exp_name"
   echo "OUTPUT (mean BV freq): ${INT_OUTFILE}"

   # Script to run
   SCRIPT_TO_RUN="BVmean_square.sh"

   # Set the brunt vaisala job parameters
   JOUT=stdout_bv_%J
   JERR=stderr_bv_%J
   JNAME="MEAN_BV"
   JQUEUE="s_medium"

   # Submit the brunt vaisala job
I   echo "Submitting the job: sh ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $nontidal_exp_path $nontidal_exp_name $sq_bvf_name $bvf_name ${INT_OUTFILE}"

   bsub -K -J $JNAME -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE \
   "sh ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $nontidal_exp_path $nontidal_exp_name $sq_bvf_name $bvf_name ${INT_OUTFILE}"

   echo "Done!"

   # Cp the outputs to the archive
   cp $WORK_DIR/${INT_OUTFILE} $OUT_DIR/

fi
if [[ $BRUNTV_FLAG == 2 ]] || [[ $BRUNTV_FLAG == 'all' ]]; then   
   # ----------------------------------
   # 2) Extract the bottom value and/or plot the output field
   echo "INPUT (mean BV freq): ${INT_OUTFILE}"
   echo "OUTPUT (bottom BV freq): ${bnbot_outfile}"

   # Script to run
   SCRIPT_TO_RUN="BVbottom_plot.py"

   # Set the bottom brunt vaisala job parameters
   JOUT=stdout_bot_%J
   JERR=stderr_bot_%J
   JNAME="BV_BOT"
   JQUEUE="s_medium"

   # Submit the bottom brunt vaisala job
   echo "Submitting the job: python ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $eas_bathy_meter $eas_mesh_mask ${INT_OUTFILE} $bvf_name ${bnbot_outfile} $bottom_bvf_name $eas_mbathy_name $eas_tmask_name $eas_lat_name $eas_lon_name $eas_bathymetry_name $order_bn $napp_bn $scheme_bn $interm_outfile $out_hrms_name $out_kbar_name"

   bsub -K -J $JNAME -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE \
   "python ${REP_SUBDIR}/${SCRIPT_TO_RUN} $WORK_DIR $eas_bathy_meter $eas_mesh_mask ${INT_OUTFILE} $bvf_name ${bnbot_outfile} $bottom_bvf_name $eas_mbathy_name $eas_tmask_name $eas_lat_name $eas_lon_name $eas_bathymetry_name $order_bn $napp_bn $scheme_bn $interm_outfile $out_hrms_name $out_kbar_name"

   echo "Done!"

   # Cp the outputs to the archive
   cp $WORK_DIR/${bnbot_outfile} $OUT_DIR/

   
   # Clean the workdir
   #rm $WORK_DIR/${INT_OUTFILE}

fi
################################### SEMIDIURNAL/DIURNAL TIDAL AMPLITUDE PERCENTAGE ##########################################

if [[ $DIURNAL_FLAG == 1 ]] || [[ $DIURNAL_FLAG == 'all' ]]; then
   
   # 1) build the file from the output of areal armonic analysis :   
   PERC_INFILE="$ampha_inpath/$ampha_infile"
   echo "INPUT  (Amplitude and phase of each tidal component): $PERC_INFILE"
   echo "OUTPUT (Semidiurnal tidal amplitude percentage field): $amperc_outfile"

   SCRIPT_TO_RUN="semidiurnal_perc.sh"
   #
   # Set the splitting job parameters
   JOUT=stdout_semid_%J
   JERR=stderr_semid_%J
   JNAME="SEMID_PERC"
   JQUEUE="s_short"

   # Submit the pre proc job
   echo "Submitting the job: sh ${REP_SUBDIR2}/${SCRIPT_TO_RUN} $WORK_DIR $PERC_INFILE $amperc_outfile $amperc_field"

   bsub -K -J $JNAME -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE \
   "sh ${REP_SUBDIR2}/${SCRIPT_TO_RUN} $WORK_DIR $PERC_INFILE $amperc_outfile $amperc_field"

   echo "Done!"

   # Cp the outputs to the archive
   cp $WORK_DIR/${amperc_outfile} $OUT_DIR/

   # 2) Plot the field
   SCRIPT_TO_RUN="semidiurnal_perc_plot.py"
   #
   # Set the splitting job parameters
   JOUT=stdout_plotSD_%J
   JERR=stderr_plotSD_%J
   JNAME="PLOT_SD"
   JQUEUE="s_short"

   # Submit the pre proc job
   echo "Submitting the job: sh ${REP_SUBDIR2}/${SCRIPT_TO_RUN} $WORK_DIR $PERC_INFILE $amperc_field $eas_lon_name $eas_lat_name $eas_bathy_meter $eas_bathymetry_name $eas_mesh_mask $eas_tmask_name"

   bsub -K -J $JNAME -P $projectid -o $EXEC_DIR/$JOUT -e $EXEC_DIR/$JERR -q $JQUEUE \
   "sh ${REP_SUBDIR2}/${SCRIPT_TO_RUN} $WORK_DIR $PERC_INFILE $amperc_field $eas_lon_name $eas_lat_name $eas_bathy_meter $eas_bathymetry_name $eas_mesh_mask $eas_tmask_name"

   echo "Done!"

fi

#############################################################################
# Check the output when possible

###############################
