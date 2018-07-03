#!/bin/bash

# This script takes an experimental design table and uses the R script "MultiPlateLayoutBlockRandomised" to generate a "pseudo-block-randomised" plate layout for it.
# This bash script is a convenience wrapper for the R script
# Please read the README file for a description of the input file

###############################################################

#set default arguments
usage="
PlateLayoutRandomistation.sh -d <designSheet> -o <outputFile> -b <batchColumnHeaders> -r <NumberOfRuns> -GH 

Options:
    -d - <string>  - designSheet - Path to experimental design file - required
    -o - <string>  - outputFileName - Output Filename prefix - optional
    -b - <string>  - batchColumns - Columns to use to create additional plots for checking batch distributions - optional
    -w - <integer> - maxWells - Maximum number of wells to use on the plate [default=96] - optional
    -G - <FLAG>    - noGenomicsControls - do not reserve wells for Genomics controls - this generally only used when the lab is doing the library prep rather than Genomics [default=FALSE] - optional
    -H - Show this help message and exit

"

#get arguments
while getopts d:o:b:w:GH opt; do
    case "$opt" in
        d) designSheet="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$designSheet" ]] ; then 
    echo "Missing/Incorrect required arguments"
    echo "$usage"
    echo "Provided arguments:"
    echo "        PlateLayoutRandomistation "$*
    exit 1
fi

# locate R script

source="${BASH_SOURCE[0]}"
source=`readlink -f ${source}`
ScriptDIR=`dirname ${source}`
ScriptFile=${ScriptDIR/bin/R}/MultiPlateLayoutBlockRandomised.R
#ScriptFile=~/Scripts/PlateLayout/R/MultiPlateLayoutBlockRandomised.R
#echo $ScriptFile
#exit

# check R installed
hash Rscript 2>/dev/null || { echo >&2 "I require R but it doesn't appear to be not installed, or is not on the PATH.  Aborting."; exit 1; }

# Run script
argumentToPass=$*
cmd="Rscript $ScriptFile $argumentToPass"


echo "Running Rscript with: $cmd"
eval $cmd
