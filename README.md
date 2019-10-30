# PlateLayout

## Installation

You can install the package from github with:

``` r
# install.packages("devtools")
devtools::install_github("crukci-bioinformatics/PlateLayout")
```

## Introduction

This package takes an experimental design table and generates a
"pseudo-block-randomised" plate layout for it. 

The package optimises the layout of the plate, and then randomises the columns
and rows. The randomisation process is carried out multiple times and the
distance between replicate samples is used to score each randomistation to see
how well distributed/badly clumpled each sample group is. The "best"
randomistation is then delivered.

The script uses one column, defined by the `primaryGroup` argument (by default
"SampleGroup"), for the optimisation/randomistaion process. If there are known
batches or other factors of interest in the experiment (e.g. different
extraction dates, sex etc.) the optimisation part of script is agnostic of
these, however, they can be considered in the scoring of each random layout via
the `batchColumns` argument.

## Inputs

The input table should be a tab separated file with samples in rows and
experimental factors and any other metadata in columns.

## Outputs

The outputs are one image (png) file for each experimental factor considered in
the layout and a table with the layout.

## Usage

### R functions

##### Randomisation

The main randomisation command is `randomizeSinglePlate`. The required input is
the metadata table as a data.frame or tibble.

By default it uses the "SampleGroup" column to optimise the layout,
however, any other column can be used instead by specifying the argument
`primaryGroup`, and in fact there are no required column headers.

Additional columns can be passed in the `batchColumns` argument, these will be
used to select the best layout from `nIter` random layouts.

The default return value is a randomised table. You can then use `plotPlate` to
plot the distribution of variables contained in the columns of the metadata 
table across the plate. e.g.:

```
library(tidyverse)
library(PlateLayout)

designSheet <- system.file("extdata", 
                           "metadata_5x3.tsv",
                           package = "PlateLayout")
bColumns <- c("ExtractionInformation", "PassageNumber") 

plateLayout <- read_tsv(designSheet) %>%  
               mutate(Replicate = str_extract(SampleName, "[0-9]$")) %>%  
               mutate(SampleGroup = str_remove(SampleName, "_[0-9]$")) %>%  
               randomizeSinglePlate(batchColumns = bColumns, nIter = 100) 

pdf("Plate_Layout.5x3.pdf")
bColumns %>%  
    map(~plotPlate(plateLayout, plotCol = .x))
dev.off()

plateLayout %>%  
    select(-SampleGroup, -RowID) %>%  
    write_tsv("Plate_Layout.5x3.tsv")
```

Alternatively, you can pass `randomizeSinglePlate` an output file prefix and it
will write out the randomized plate layout as a tsv file and individual plots
for each of the `batchColumns`.

```
designSheet <- system.file("extdata", 
                           "metadata_12x3.tsv",
                           package = "PlateLayout")
bColumns <- c("ExtractionInformation", "PassageNumber") 
outputFile <- "Plate_Layout.12x3"

read_tsv(designSheet) %>%  
    randomizeSinglePlate(outputFile = outputFile, 
                         batchColumns = bColumns, 
                         nIter = 100) 
```

**Note:** by default the script runs in parallel across 4 cores. You can change
this with the `nCores` argument.

#### CRUK CI Bioinformatics metadata template

There is a utility script `fixMetaForm` that will load a CRUK CI Bioinformatics
Core metadata template Excel file and tidy it into a table suitable for use
with the above randomisation function.

### Bash script

In the *bin* directory of the Git repository is a bash script
`RandomizeSinglePlate.sh` that can be used at the command line to run the
randomisation. e.g.

```
RandomizeSinglePlate.sh -d example_metadata.tsv \
                        -o Plate_Layout \
                        -b SampleType,Patient \
                        -r 100
```

Please use `RandomizeSinglePlate.sh -H` for usage.

## Plate Layout Process

The algorithm for randomising the plate layout is as follows:

### Optimisation
1. Randomise the order of the sample groups and the order of the replicates
   within each sample group:  
    a) assign each sample group a random number from 1-nGroups  
    b) assign each sample a random number from 1-nSamples  
    c) sort the table by group and then sample  

2. Assign each sample to a well by traversing the plate diagonally:    
    a) determine the number of columns that will be needed to accomodate all the
samples on the plate    
    b) the first sample in the table is then assigned to A1, the second to B2
      and so on. On reaching the final column continue the next row from column
      1, and continue the next column from Row A when Row H has been reached.  
    The purpose of assigning wells in this way is to distribute the
      replicates from each sample group across the rows and the columns so
      that, as much as possible, replicates are not in the same row or column as
      others from the same sample group.  

### Randomisation

5. Randomise rows and columns as blocks - this maintains the distribution of
   replicates relative to each other, but randomises the overall plate layout.  
    a) each row is assigned a randomn number from 1-8 and the rows are re-ordered
     according to this  
    b) each column is assigned a randomn number from 1-nColumns and the columns
      are ordered according to this  

    e.g. for 5 rows with 40 samples:
    
        Step 4 ->
             1  17  33   9  25
            26   2  18  34  10
            11  27   3  19  35
            36  12  28   4  20
            21  37  13  29   5
             6  22  38  14  30
            31   7  23  39  15
            16  32   8  24  40

        Step 5i ->
            31   7  23  39  15
            11  27   3  19  35
            21  37  13  29   5
            16  32   8  24  40
            36  12  28   4  20
             1  17  33   9  25
            26   2  18  34  10
             6  22  38  14  30

        Step 5ii ->
            39  31  23  15   7
            19  11   3  35  27
            29  21  13   5  37
            24  16   8  40  32
             4  36  28  20  12
             9   1  33  25  17
            34  26  18  10   2
            14   6  38  30  22
            
6. Score the randomised plate for distribution of each factor:  
    A score for each sample is calculated as the sum of 1 / (Euclidean Distance
    to other samples of the same type where dist is <= 2).  
    The distribution score for the plate is the sum of scores for all samples
    for all "batch columns".
                    
7. The randomisation process is carried out 10000 times and the layout with the
   lowest "distribution score" is accepted.
