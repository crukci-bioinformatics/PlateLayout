#' Randomize Single Plate
#'
#' Randomize up to 96 samples across a single 96 well plate. 
#' 
#' @param designTable The meta data table - a data.frame/tibble.
#' @param outputFile The root name for the output files if required.
#' @param primaryGroup The name of a column in the metadata table to use for
#' optimising the distribution of samples across the plate
#' @param batchColumns A character vector wit the names of additional columns
#' in the metadata table to use in assessing the distribution of the samples
#' across the plate.
#' @param nIter The number of random layouts to generate and assess.
#' @param nCores The number of cores to use. Set to 1 to run serial.
#' @details The metadata (\code{designTable}) should be a data.frame/tibble. The 
#' column names can be anything. By default the function expects one column 
#' named \code{SampleGroup} that it will use to distribute the samples across
#' the plate. This can be set to any other column using the \code{primaryGroup}
#' argument.
#' @details The columns specified in \code{batchColumns} are used in addition
#' to \code{primaryGroup} to assess which of the \code{nIter} layouts is 
#' optimal in terms of the distribution of samples. This is assessed based on
#' the sum of Euclidean distances between replicate samples within each column.
#' @return If no \code{outputFile} is provided, a tibble containing the original
#' metadata and the well assignments for each sample. If an \code{outputFile} is
#' provided, nothing - the function writes out the plate layout to a tsv file 
#' and a plot for each of `primaryGroup` and `batchColumns` to png files.
#' @examples
#' library(readr)
#' designSheet <- system.file("extdata", "metadata_12x3.tsv",
#'                            package = "PlateLayout")
#' bColumns <- c("ExtractionInformation", "PassageNumber") 
#' designTable <- read_tsv(designSheet)
#' # Create an R object in the session
#' plateLayout <- randomizeSinglePlate(designTable,  batchColumns = bColumns,
#'                      nIter = 100) 
#'
#' # Or directly output a table and plots
#' outputFile <- "Test"
#' randomizeSinglePlate(designTable, 
#'                      outputFile = outputFile, 
#'                      batchColumns = bColumns,
#'                      nIter = 100) 
#' @export
randomizeSinglePlate <- function(designTable, 
                                 outputFile =  NULL, 
                                 primaryGroup = "SampleGroup", 
                                 batchColumns = NULL, 
                                 nIter = 10000,
                                 nCores = 4){
    # Set parallel processing
    options(future.supportsMulticore.unstable="quiet") # turn off warning about forking in RStudio
    if(nCores>1){ plan(multicore, workers = nCores) }

    # Checks
    if(length(primaryGroup) > 1) {
        stop("Only specify 1 column for the 'primaryGroup'")
    }
    if(!primaryGroup%in%colnames(designTable)){
        stop("'", primaryGroup, "' is not a column in the design sheet")
    }
    if(any(!batchColumns%in%colnames(designTable))){
       stop("Not all of the provided 'batchColumns' are in the design sheet")
    }
    if(nrow(designTable) > 96){
        stop("There are too many samples in the sample sheet for a single plate")
    }
    if(nrow(designTable) > 94){
        warning("There are more than 94 samples. Are you sure you do not ",
                "want to leave wells for the Genomics positive controls?")
    }

    # Modify batch columns
    bCols <- c(primaryGroup, batchColumns) %>%  
        str_subset("Replicate", negate = TRUE) %>%  
        unique()

    message("Run ", nIter, " random layouts")
    layouts <- future_map(1:nIter, 
                          ~plateRandomisation(designTable, 
                                              primaryGroup = primaryGroup),
                          .progress = TRUE)

    message("Calculate distribution scores for each layout and find best layout")
    bestLayout <- getMinScore(layouts, bCols)
        
    # Choose best layout
    finalLayout <- pluck(layouts, bestLayout)

    if(!is.null(outputFile)){
        message("Output plots and layout")
        # Save a plot for each batch column
        bCols %>%  
            walk(outputPlatePlot, datTab = finalLayout, outFileName = outputFile)    

        # Write out the layout
        finalLayout %>%  
            select(-RowID) %>%  
            write_tsv(str_c(outputFile, ".PlateLayout.tsv"))
    }else{
        return(finalLayout)
    }
}
