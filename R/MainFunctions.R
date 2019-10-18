#' Randomize Single Plate
#'
#' Randomize up to 96 samples across a single 96 well plate. 
#' 
#' @param designSheet The path to a tab delimited sample metadata file.
#' @param outputFile The root name for the various output files.
#' @param primaryGroup The name of a column in the metadata table to use for
#' optimising the distribution of samples across the plate
#' @param batchColumns A character vector wit the names of additional columns
#' in the metadata table to use in assessing the distribution of the samples
#' across the plate.
#' @param nIter The number of random layouts to generate and assess.
#' @param nCores The number of cores to use. Set to 1 to run serial.
#' @details The sample metadata file should be a tab delimited table. The 
#' column names can be anything. By default the function expects one column 
#' named \code{SampleGroup} that it will use to distribute the samples across
#' the plate. This can be set to any other column using the \code{primaryGroup}
#' argument.
#' @details The columns specified in \code{batchColumns} are used in addition
#' to \code{primaryGroup} to assess which of the \code{nIter} layouts is 
#' optimal in terms of the distribution of samples. This assessed based on the
#' sum of Euclidean distances between replicate samples within each column.
#' @return Nothing. The function writes out the plate layout to a tsv file and
#' a plot for each of `primaryGroup` and `batchColumns` to png files.
#' @examples
#' designSheet <- system.file("extdata", "metadata_template_TEST_12x3.tsv",
#'                            package = "PlateLayout")
#' outputFile <- "Test"
#' bColumns <- c("ExtractionInformation", "PassageNumber") 
#' randomizeSinglePlate(designSheet, outputFile, batchColumns = bColumns,
#'                      nIter = 100) 
#' @export
randomizeSinglePlate <- function(designSheet, 
                                 outputFile, 
                                 primaryGroup = "SampleGroup", 
                                 batchColumns = NULL, 
                                 nIter = 10000,
                                 nCores = 4){
    # Set parallel processing
    options(future.supportsMulticore.unstable="quiet") # turn off warning about forking in RStudio
    if(nCores>1){ plan(multiprocess, workers = nCores) }

    # Checks
    if(!file.exists(designSheet)){ stop("Cannot locate ", designSheet) }
    if(length(primaryGroup) > 1) {
        stop("Only specify 1 column for the 'primaryGroup'")
    }

    # read sample sheet
    samsht <- read_tsv(designSheet, col_types = cols(.default = "c"))
    if(!primaryGroup%in%colnames(samsht)){
        stop("'", primaryGroup, "' is not a column in the design sheet")
    }
    if(any(!batchColumns%in%colnames(samsht))){
       stop("Not all of the provided 'batchColumns' are in the design sheet")
    }
    if(nrow(samsht) > 96){
        stop("There are too many samples in the sample sheet for a single plate")
    }
    if(nrow(samsht) > 94){
        warning("There are more than 94 samples. Are you sure you do not ",
                "want to leave wells for the Genomics positive controls?"):x
    }

    # Modify batch columns
    bCols <- c(primaryGroup, batchColumns) %>%  
        str_subset("Replicate", negate = TRUE) %>%  
        unique()

    message("Run ", nIter, " random layouts")
    layouts <- future_map(1:nIter, 
                          ~plateRandomisation(samsht, 
                                              primaryGroup = primaryGroup),
                          .progress = TRUE)

    message("Calculate distribution scores for each layout and find best layout")
    bestLayout <- getMinScore(layouts, bCols)
        
    # Choose best layout
    finalLayout <- pluck(layouts, bestLayout)

    message("Output plots and layout")
    # Save a plot for each batch column
    bCols %>%  
        walk(savePlots, datTab = finalLayout, outFileName = outputFile)    

    # Write out the layout
    finalLayout %>%  
        select(-RowID) %>%  
        write_tsv(str_c(outputFile, ".PlateLayout.tsv"))
}
