library(optparse)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))

#Input options
options_list <- list(
    make_option(c('--designSheet', '-d'), type='character', help="Path to
                experimental design file - required", dest="designSheet",
                default="<undefined>"),
    make_option(c('--output', '-o'), type='character', help="Output Filename -
                optional", dest="outputFile", default="<undefined>"),
    make_option(c('--batchColumns', '-b'), type='character', help="Columns to
                be used to assess distribution of samples across plate.",
                dest="batchColumns", default=NULL),
    make_option(c('--nIter', '-n'), type='integer', help="The number of random
                layouts to test.", dest="nIter", default=10000),
    make_option(c('--nCores', '-c'), type='integer', help="The number of cores
                to use in parallel.", dest="nCores", default=4)
                )
parser <- OptionParser(option_list=options_list, usage="%prog [options]")
opts <- parse_args(parser)
designFile <- opts$designSheet
outputFile <- opts$outputFile
batchColumns <- opts$batchColumns 
nIter <- opts$nIter 
nCores <- opts$nCores 

designSheet <- "../test/metadata_template_TEST_12x3.csv"
outputFile <- "Test"
batchColumns <-"ExtractionInformation,PassageNumber" 
nIter <- 100
nCores <- 4

# Set parallel processing
plan(multiprocess, workers = nCores)

# Load functions
source("RandomisePlate.R")
source("utilities.R")

# read sample sheet
samsht <- read_csv(designSheet, col_types = cols(.default = "c")) 

# Modify batch columns
bCols <- batchColumns %>%  
    str_split(",", simplify = TRUE) %>%  
    c("SampleGroup") %>%  
    str_subset("Replicate", negate = TRUE) %>%  
    unique()

message("Run ", nIter, " random layouts")
layouts <- future_map(1:nIter, ~runRandomisation(samsht), .progress = TRUE)

message("Calculate distribution scores for each layout and find best layout")
bestLayout <- getMinScore(layouts, bCols)
    
# Choose best layout
finalLayout <- pluck(layouts, bestLayout)

message("Output plots and layout")
# Save a plot for each batch column
bCols %>%  
    walk(savePlots, datTab = finalLayout, outFileName = outputFile)    

# Write out the layout

write_tsv(finalLayout, str_c(outputFile, ".PlateLayout.tsv"))

