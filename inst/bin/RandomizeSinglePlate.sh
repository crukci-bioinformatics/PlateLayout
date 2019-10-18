#!/bin/bash

# This script takes an experimental design table and uses the R package
# PlateLayout to generate a "pseudo-block-randomised" plate layout for it.
# Please read the package help for a description of the input file

###############################################################

#set default arguments
usage="
PlateLayoutRandomistation.sh -d <designSheet> -o <outputFile> -b
<batchColumnHeaders> -r <NumberOfRuns> -GH 

Options:
    -d - <string>  - designSheet - Path to experimental design file - required
    -o - <string>  - outputFileName - Output Filename prefix
    -p - <string>  - primaryGroup - Column to use for optimisation. Default
                     'SampleGroup' - optional
    -b - <string>  - batchColumns - Columns to use to create additional plots
                     for checking batch distributions. Separate multiple
                     columns with commas. No kneed to specify 'SampleGroup'
                     - optional
    -r - <integer> - Number of iterations to run. Default 10000 - optional
    -n - <integer> - Number of cores to use. Default 4. - optional
    -H - Show this help message and exit

"

#get arguments
while getopts d:o:p:b:r:n:GH opt; do
    case "$opt" in
        d) designSheet="$OPTARG";;
        o) outputPrefix="$OPTARG";;
        p) primaryGroup="$OPTARG";;
        b) batchColumns="$OPTARG";;
        r) nIter="$OPTARG";;
        n) nCores="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$designSheet" ]] || [[ ! ${outputPrefix} ]] ; then 
    echo "Missing/Incorrect required arguments"
    echo "$usage"
    echo "Provided arguments:"
    echo "        PlateLayoutRandomistation "$*
    exit 1
fi

# check R installed
hash Rscript 2>/dev/null || \
        { echo >&2 "R is not on the PATH.  Aborting."; exit 1; }

# Run script
cmd="Rscript -e 'library(PlateLayout); randomizeSinglePlate(\"${designSheet}\", 
    outputFile=\"${outputPrefix}\""
if [[ ${primaryGroup} ]]; then cmd=${cmd}", primaryGroup=\"${primaryGroup}\""; fi
if [[ ${batchColumns} ]]; then 
    bColVector=`echo ${batchColumns} | sed -e 's/^/c(\"/' -e 's/,/\", \"/g' -e 's/$/\")/'`
    cmd=${cmd}", batchColumns=${bColVector}"
fi
if [[ ${nIter} ]]; then cmd=${cmd}", nIter=${nIter}"; fi
if [[ ${nCores} ]]; then cmd=${cmd}", nCores=${nCores}"; fi
cmd=${cmd}"); sessionInfo()'"

echo "Running Rscript with: "
echo "    $cmd"
eval $cmd

