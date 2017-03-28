# PlateLayout

Plate randomistation using R script MultiPlateLayoutBlockRandomised  
Use the shell script in the `bin` directory to access the script

----

usage:  
`PlateLayoutRandomisation -d <designSheet> -o <outputFile> -b <batchColumnHeaders> -w <maxWells> -G -H`

Options:  
>    -d - <*string*>  - designSheet - Path to experimental design file - **required**  
>    -o - <*string*>  - outputFileName - Output Filename prefix - **optional**  
>    -b - <*string*>  - batchColumns - Columns to use to create additional plots for checking batch distributions - **optional**  
>    -w - <*integer*> - maxWells - Maximum number of wells to use on the plate [default=96] - **optional**  
>    -G - <*FLAG*>    - noGenomicsControls - do not reserve wells for Genomics controls - this is generally only used when the lab is doing the library prep rather than Genomics [default=FALSE] - **optional**  
>    -H - <*FLAG*> - Show help message and exit

## INTRODUCTION

This script takes an experimental design table and generates a "pseudo-block-randomised" plate layout for it. 

Whenever the number of samples on a plate is not divisable by 8, the script will add up to two Genomic Controls and as many water controls as necessary to ensure that columns on the plate are complete.

The script first optimises the layout of the plate, and then randomises the columns and rows. The randomisation process is carried out 10000 times and an algorithm is used to score each randomistation to see how well distributed/badly clumpled each sample group is. The "best" randomistation is the delivered.

If there are more than 94 samples the script will generate layouts for multiple plates to accomodate all the samples - the last two wells on the 96 well plate are reserved for genomic controls. Samples will be randomly spread between the plates, ensuring even distribution of different sample groups to avoid batch effects. A pdf will be created with bar plots showing the distributions of Sample groups and replicates between the plates.

The script requires columns called "SampleGroup" and "Replicate" that will be used for the optimisation/randomistaion process. If there are known batches or other factors of interest in the experiment (e.g. different extraction dates, sex etc.) the optimisation part of script is agnostic of these, but the ranodmisation can take account of them. This is achieved by providing the `-b` option with a comma separated string of column headers for which plots are required e.g. `-b ExtractDate,Sex,PassageNumber`.

The script generates layout plots for each factor and, in the case of multi-plate experiments, additional bar plots showing the distribution of these characteristics across the plates.  These plots can then be used to assess if the characteristics are evenly distributed across and between the plates or if the script needs to be re-run. 

For some library preps sample pool size is limited to less than 96 (e.g. methylation is currently 24). In this case the `-w` option allows you to reduce the maximum number of wells that can be used on a single plate.

The `-G` will cause the script to not reserve wells for the genomics controls. These controls are added by Genomics in order to verify the library prep in cases where the sequecing has failed due to problems caused in sample preparation/extraction. 

## REQUIRMENTS

An installation of R is required:  
  [CRAN](https://cran.r-project.org/)
  
In addition, the following R packages are required:
  * dplyr
  * ggplot2
  * RColorBrewer
  * optparse

## INPUTS

The input table should be a csv file with a minimum of three columns headed "SampleName", "SampleGroup" and "Replicate". Any additional columns will be ignored unless specified with `-b`.

## OUTPUTS

The outputs are at least one image file for each plate showing the layout, indicating both sample group and replicate number, a table with the layout, and, in multi-plate experiments, a pdf with multiple barplots showing the distribution of different sample characteristics between the plates.

The output table will include all the data in the input table, plus 3 or 4 additional columns indicating the plate number (if multiple plates), the row of the well (A-H), the column of the well (1-12) and the coordinate of the well (A1, A2...H7,H8).

## PLATE LAYOUT PROCESS

The algorithm for randomising the plate layout is as follows:

### OPTIMISATION
1. Randomise the order of the sample groups and the order of the replicates within each sample group:
    a) assign each sample group a random number from 1-nGroups
    b) assign each sample a random number from 1-nSamples
    c) sort the table by group and then sample

2. Assign each sample to a plate

 **for each plate...**

3. Add any genomics controls and water controls. Genomics controls are added to end of the table. Water controls are added randomly within the table.

4. Assign each sample to a well by traversing the plate diagonally:  
   i. determine the number of columns that will be needed to accomodate all the samples on the plate  
  ii. the first sample in the table is then assigned to A1, the second to B2 and so on. On reaching the final column continue the next row from column 1, and continue the next column from Row A when Row H has been reached.  
    * If a well has already been used, shift over 1 column - this is because the method hits problems if there is an even number of columns as the pattern cycles (e.g. with 8 columns, it would just keep using A1, B2, C3, D4, E5, F6, G7, H8)  
    * The purpose of assigning wells in this way is to distribute the replicates from each sample group across the rows and the columns so that, as much as possible, replicates are not in the same row or column as others from the same sample group.  

### Randomisation

5. Randomise rows and columns as blocks - this maintains the distribution of replicates relative to each other, but randomises the overall plate layout.
  i. each row is assigned a randomn number from 1-8 and the rows are re-ordered according to this
  ii. each column is assigned a randomn number from 1-nColumns and the columns are ordered according to this
  iii. Genomics controls are switched with the samples in the last 1-2 wells (e.g. for 5 columns with wells E7 and E8 if there are two genomics controls)

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
    "Sistribution Score" = Sum of:
                        For each well:
                            For each factor:
                                Score for nearby samples of the sample type:
                                    Score 2 for each sample in the four adjacent wells that of the sample group
                                    Score 1 for each sample in the six next nearest wells
                                    
                                    |---|---|---|---|---|
                                    |   |   | 1 |   |   |
                                    |---|---|---|---|---|
                                    |   | 1 | 2 | 1 |   |
                                    |---|---|---|---|---|
                                    | 1 | 2 | X | 2 | 1 |
                                    |---|---|---|---|---|
                                    |   | 1 | 2 | 1 |   |
                                    |---|---|---|---|---|
                                    |   |   | 1 |   |   |
                                    |---|---|---|---|---|
                    
7. The randomisation process is carried out 10000 times and the layout with the lowest "distribution score" is accepted.

TEST DATA
=========

The folder test contains 11 example design tables, which can be use to assess the behaviour of the script. The script file `TestScript.sh` contains some examples of how to use the various options for the script. Run `TestScript.sh` in a local directory to get example outputs.
