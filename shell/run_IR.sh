#!/bin/bash

# Script runs index_repeat from command line
# See relevant documentation in function (/project/souce/dev/...)
# Output filename is always placed in: /project/analysis/outputs/IR

function usage()
{
 echo "Usage: run_IR.sh [ infile ] [ step ] [ no. grains ] [ repeat ] [ seed ] [ crystal ] [ index ] [ bin size (md only) ] [ binning type (md only) [ output name ]" 
 echo "Usage: alternatively, will run interactively if no arguments given"
}


#------------------------------------------------------------------------
# check if running quietly or accepting inputs

if [[ $# -eq 0 ]] 
then 

  echo
  printf "Input file:........ "
  read infile
  printf "Step:.............. "
  read step
  printf "No. grains:........ "
  read n
  printf "Repeat:............ "
  read repeat
  printf "Seed:.............. "
  read seed
  printf "Crystal............ "
  read crystal
  printf "Index:............. "
  read index

  # if we're running m_indexDisc, need to know bin size
  case $index in
     md|MD) printf "Bin size (deg)..... "
            read bin
            printf "Binning type....... "
            read binning
            
             needBin=1 # we do have bin size
             ;;
          *) needBin=0 # we don't need bin size
             ;;
  esac

  printf "Output file name:.. "
  read outname
  echo "Running function with user inputs..."
  echo

else

  # check if discrete m index (so we need bin size)
  case $7 in
     md|MD)  [[ $# -ne 10 ]] && usage && exit 1
             infile=$1
             step=$2
             n=$3
             repeat=$4
             seed=$5
             crystal=$6
             index=$7
             bin=$8
             binning=$9
             outname=${10}
             needBin=1 # we do need bin 
             ;;

  # index must be either j or m cont.
          *) [[ $# -ne 8 ]] && usage && exit 1
             infile=$1
             step=$2
             n=$3
             repeat=$4
             seed=$5
             crystal=$6
             index=$7
             outname=$8
             needBin=0 # we dont need bin
             ;;
  esac

fi

#------------------------------------------------------------------------
# setup important dirs

#### CHANGE TO RELEVANT DIRECTORIES ####
outdir="/nfs/see-fs-01_teaching/ee12lmb/project/analysis/outputs/IR"
devdir="/nfs/see-fs-01_teaching/ee12lmb/project/source/dev"
outfile=$outdir/$outname
#echo "OUTNAME: $outname"
#echo "OUTPATH: $outfile"

#------------------------------------------------------------------------
# input checks
[[ ! -f $infile ]] && echo "Input file not found!" && usage && exit 1
[[ -f $outfile ]] && echo "Output file already exists!" && usage && exit 1

#------------------------------------------------------------------------
# Run matlab function


##### CHANGE PATH IN ADDPATH COMMAND TO THE TEXTURE STRENGTH PACKAGE LOACTION IN BOTH COMMANDS BELOW #####
if [[ $needBin -eq 0 ]]
then
  matlab -nodesktop -nodisplay -nosplash -r "addpath('/nfs/see-fs-01_teaching/ee12lmb/project/source/dev/'); setup_env; index_repeat('$infile',$step,$n,$repeat,$seed,'crystal','$crystal','index','$index','outfile','$outfile','-v'); exit;"

elif [[ $needBin -eq 1 ]]
then 
  matlab -nodesktop -nodisplay -nosplash -r "addpath('/nfs/see-fs-01_teaching/ee12lmb/project/source/dev/'); setup_env; index_repeat('$infile',$step,$n,$repeat,$seed,'crystal','$crystal','index','$index','bin',$bin,'binning','$binning','outfile','$outfile','-v'); exit;"
fi

exit 0
