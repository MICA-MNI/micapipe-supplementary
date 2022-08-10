#!/bin/bash
# MICA T1w/T2
#
# This script will create a T1w/T2 using two files from the rawdata directory
#
#
#
#
#

help() {
echo -e "
\033[38;5;141mCOMMAND:\033[0m
   $(basename $0)

\033[38;5;141mARGUMENTS:\033[0m
\t\033[38;5;197m-sub\033[0m 	        : Subject identification (no 'sub-')
\t\033[38;5;197m-out\033[0m      	: Output directory for the processed files <anat ORIG and face-warps>.
\t\033[38;5;197m-bids\033[0m      	: Path to BIDS directory
\t\033[38;5;120m-ses\033[0m 	        : OPTIONAL flag that indicates the session name (if omitted will manage as SINGLE session)

\033[38;5;141mOPTIONAL ARGUMENTS:\033[0m
\t\033[38;5;197m-t2Str\033[0m 	    : T2 image string
\t\033[38;5;197m-threads\033[0m        : Number of threads (Default is 6)

\033[38;5;141mUSAGE:\033[0m
\033[38;5;141m$(basename $0)\033[0m  \033[38;5;197m-sub\033[0m <subject_id> \033[38;5;197m-out\033[0m <outputDirectory> \033[38;5;197m-bids\033[0m <BIDS-directory>\n \033[38;5;197m-all\033[0m

\033[38;5;141mDEPENDENCIES:\033[0m
    > ANTs        2.3.3   (https://github.com/ANTsX/ANTs)
    > FSL         6.0     (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
    > AFNI        20.3.03 (https://afni.nimh.nih.gov/download)

RRC
McGill University, MNI, MICA-lab, October 2021
https://github.com/MICA-MNI
http://mica-mni.github.io/
"
}

# source print funtions for MICAPIPE
source "${MICAPIPE}/functions/utilities.sh"

# -----------------------------------------------------------------------------------------------#
#			ARGUMENTS
# Create VARIABLES
for arg in "$@"
do
  case "$arg" in
  -h|-help)
    help
    exit 1
  ;;
  -sub)
    id=$2
    shift;shift
  ;;
  -out)
    out=$2
    shift;shift
  ;;
  -bids)
    BIDS=$2
    shift;shift
  ;;
  -ses)
    SES=$2
    shift;shift
  ;;
  -t2Str)
    t2wStr=$2
    shift;shift
  ;;
  -threads)
    threads=$2
    shift;shift
  ;;
  -*)
    Error "Unknown option ${2}"
    help
    exit
  ;;
    esac
done

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# argument check out & WARNINGS
arg=($id $out $BIDS)
if [ "${#arg[@]}" -lt 3 ]; then
Error "One or more mandatory arguments are missing:
               -sub  : $id
               -out  : $out
               -bids : $BIDS"
help; exit 1; fi

# Get the real path of the Inputs
out=$(realpath $out)/micapipe
BIDS=$(realpath $BIDS)
id=${id/sub-/}

# Number of session (Default is "ses-pre")
if [ -z ${SES} ]; then SES="ses-pre"; else SES="ses-${SES/ses-/}"; fi

# THREADS should be defined as GLOBAL variable maybe
if [[ -z $threads ]]; then export threads=10; fi

# Assigns BIDS variables names
bids_variables $BIDS $id $out $SES

# Check BIDS directory
if [ ! -d ${subject_bids} ]; then Error "$id was not found on the BIDS directory\n\t     ${subject_bids}"; exit 0; fi

# Check inputs: Nativepro T1
T1nativepro="${proc_struct}/${idBIDS}"_space-nativepro_t1w.nii.gz
if [ ! -f "${T1nativepro}" ]; then Error "Subject $id doesn't have T1_nativepro: ${T1nativepro}"; exit; fi

# T1w image to process (DEFAULT is *T2w.nii*)
if [ -z ${t2wStr} ]; then t2wStr=DEFAULT; else t2wStr="$t2wStr"; fi
Note "T2         =" "$t2wStr"
# Manage manual inputs: T1w images
if [[ "$t2wStr" != "DEFAULT" ]]; then
  IFS=',' read -ra bids_t2wStr <<< "$t2wStr"
  for i in "${!bids_t2wStr[@]}"; do bids_t2wStr[i]=$(ls "${subject_bids}/anat/${idBIDS}_${bids_t2wStr[$i]}.nii"* 2>/dev/null); done
  bids_T2ws=("${bids_t2wStr[@]}")
  N="${#bids_t2wStr[*]}"
  Info "Manually selected T2 weighted images: $t2wStr, N=${N}"
else
  bids_T2ws=($(ls "$subject_bids"/anat/*T2*nii 2>/dev/null))
fi

Title "Running MICA T1w/T2"
Note "Subject         =" "$id"
Note "Session         =" "${SES/ses-/}"
Note "T2w         =" "${bids_T2ws[0]}"
# Check if T2 file exist
if [ ! -f "${bids_T2ws[0]}" ]; then Error "Subject $id doesn't have an associated T2 weighted image: ${bids_T2ws[0]}"; exit; fi

# -----------------------------------------------------------------------------------------------#
#			  Timer & Beginning
aloita=$(date +%s)

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporal directory
tmp=${tmp}/${RANDOM}_micapipe_t1t2_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

here=$(pwd)
#------------------------------------------------------------------------------#
Note "T1w         =" "$T1nativepro"
Info "micapipe will use $threads threads for multicore processing"
Info "Saving temporal dir: $nocleanup"
Info "ANTs will use $threads threads"

# Creates OUT directory if it doesn't exist
Info "Output directory: $out"


# Variables
t1wt2w="${proc_struct}/${idBIDS}_space-nativepro_t1wt2w.nii.gz"
t2w_tmp="${tmp}/${idBIDS}_space-nativepro_t2w.nii.gz"
t2w_res="${tmp}/${idBIDS}_space-nativepro_t2w_rescaled.nii.gz"

# -----------------------------------------------------------------------------------------------
Info "Affine registration between T2 and T1_nativepro"
T2w_2_natT1w="${dir_warp}/${idBIDS}"_from-T2w_to-nativepro_mode-image_desc-affine
T2w_2_natT1w_mat=${T2w_2_natT1w}_0GenericAffine.mat

Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro -m ${bids_T2ws[0]} -o ${T2w_2_natT1w}_ -t a -n $threads -p f
if [ ! -f ${T2w_2_natT1w_mat} ]; then Error "Affine registration failed"; cd $here;  Do_cmd exit 0; fi

# Apply the warpflied
Do_cmd antsApplyTransforms -d 3 -u int -i ${bids_T2ws[0]} -r $T1nativepro -n NearestNeighbor \
        -t $T2w_2_natT1w_mat -o $t2w_tmp

Info "Generate the T1w over the T2 weighted image"
# Intensity Non-uniform correction - N4
Do_cmd N4BiasFieldCorrection  -d 3 -i "$t2w_tmp" -r -o "$t2w_tmp" -v

# Rescale intensity [100,0]
Do_cmd ImageMath 3 "$t2w_res" RescaleImage "$t2w_tmp" 0 100

# Generate the T1/T2w
Do_cmd ImageMath 3 "$t1wt2w" / "$T1nativepro" "$t2w_res"

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
Do_cmd cd $here
if [[ -z $nocleanup ]]; then Do_cmd rm -rf $tmp ${out_warps}/anat; else Info "tmp directory was not erased: ${tmp}"; fi
Note "Temporal directory:  " $tmp
Note "Output directory:    " $out

# -----------------------------------------------------------------------------------------------#
#			 Total Running Time
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print $eri/60 | perl)
Title "TOTAL running time:\033[38;5;220m $(printf "%0.3f\n" ${eri}) minutes \033[38;5;141m"
bids_variables_unset
