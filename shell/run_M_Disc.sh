#!/bin/bash

# Script runs m_indexDisc from command line
# See relevant documentation in function (/project/souce/dev/...)
# Output filename is always placed in: /project/analysis/outputs/M_Disc

function usage()
{
 echo "Usage: run_M_Disc.sh [ infile ] [ no. grains ] [ seed ] [ crystal ] [ bin ] [ binning ] [ output name ]" 
 echo "Usage: alternatively, will run interactively if no arguments given"
}


#------------------------------------------------------------------------
# check if running quietly or accepting inputs

if [[ $# -eq 0 ]] 
then 

  echo
  printf "Input file:........ "
  read infile
  printf "No. grains:........ "
  read n
  printf "Seed:.............. "
  read seed
  printf "Crystal............ "
  read crystal
  printf "Bin size (deg)..... "
  read bin
  printf "Binning type....... "
  read binning
  printf "Output file name:.. "
  read outname
  echo "Running function with user inputs..."
  echo

else

  [[ $# -ne 7 ]] && usage && exit 1
  infile=$1
  n=$2
  seed=$3
  crystal=$4
  bin=$5
  binning=$6
  outname=$7

fi

#------------------------------------------------------------------------
# setup important dirs

#### CHANGE TO RELEVANT DIRECTORIES ####
outdir="/nfs/see-fs-01_teaching/ee12lmb/project/analysis/outputs/M_Disc"
devdir="/nfs/see-fs-01_teaching/ee12lmb/project/source/dev/"
outfile=$outdir/$outname

#------------------------------------------------------------------------
# input checks
[[ ! -f $infile ]] && echo "Input file not found!" && exit 1
[[ -f $outfile ]] && echo "Output file $outfile already exists!" && exit 1

#------------------------------------------------------------------------
# Run matlab function

##### CHANGE PATH IN ADDPATH COMMOAND BELOW TO TEXTURE STRENGTH PACKAGE LOCATION #####
matlab -nodesktop -nodisplay -nosplash -r "addpath('$devdir'); setup_env; m_indexDisc('$infile',$n,$seed,'crystal','$crystal','bin',$bin,'binning','$binning','outfile','$outfile'); exit;"

exit 0
