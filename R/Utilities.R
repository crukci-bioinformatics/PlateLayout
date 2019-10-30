# Test the distribution of a samples in a single column across the plate
# Use the sum of 1/(Euclidean distances between samples less than or equal to 2)
# For a given sample: 
#     0   0   0.5   0   0
#     0 0.707  1  0.707  0
#    0.5   1   X    1  0.5 
#     0 0.707  1  0.707  0
#     0   0   0.5   0   0

getScore <- function(rowCols){
    dst <- dist(rowCols, method = "euclidean")
    sum(1 / dst[dst <= 2])
}

calcDistrScore <- function(batchColumn, dat){
    dat %>%  
        select(bCol = batchColumn, Column, RowID) %>%  
        group_by(bCol) %>%  
        nest(dat = c(Column, RowID)) %>%  
        mutate(score = map_dbl(dat, getScore)) %>%  
        pull(score) %>%  
        sum()} 

# get distr scores for each batch column of a sample sheet 

getDistrScores <- function(samSht, batchColumns){ 
    batchColumns %>%  
        set_names(batchColumns) %>%  
        future_map_dfr(calcDistrScore, samSht)
}    

# get scores for all random layouts for all batch columns, normalised within
# each batch, sum, find best layout as minimum score

normScores <- function(scr){ scr / median(scr) }

getMinScore <- function(layoutList, batchCols){
    future_map_dfr(layoutList, getDistrScores, batchColumns = batchCols, 
                          .progress = TRUE) %>% 
        mutate_all(normScores) %>%  
        rowSums() %>%  
        which.min()
}

###############################################################################

# create a colour scheme using color brewer, few options depending on number of
# groups
getCols <- function(dat, plotCol){
    len <- pull(dat, plotCol) %>% unique() %>% length()
    if(len<=2){
        wcols <- c("#00B6ED", "#EC008C")
    }else if(len<=12){
        wcols <- RColorBrewer::brewer.pal(len, "Set3")
    }else if(len<=21){
        wcols <- c(RColorBrewer::brewer.pal(12, "Set3"), 
                   suppressWarnings(
                      RColorBrewer::brewer.pal(len - 12, "Set1")
                      ))
    }else{
        wcols <- sample(rainbow(len), len)
    }
    wcols
}

#' Plot Plate Layout
#'
#' Plot the plate with wells coloured according to a column in the samplesheet,
#' if there is column called Replicate the numbers will also be printed
#' 
#' @param dat The plate layout tables
#' @param plotCol The column of the table by which to colour the wells
#' @param wellCols A named vector of columns for colouring the wells. Names
#' should match the sample groups in plotCol
#' @return A ggplot object
#' @examples
#' library(readr)
#' designSheet <- system.file("extdata", "metadata_12x3.tsv",
#'                            package = "PlateLayout")
#' bColumns <- c("ExtractionInformation", "PassageNumber") 
#' designTable <- read_tsv(designSheet)
#' # Create an R object in the session
#' plateLayout <- randomizeSinglePlate(designTable,  
#'                                     batchColumns = bColumns,
#'                                     nIter = 100) 
#' plotPlate(plateLayout)
#' plotPlate(plateLayout, "ExtractionInformation")
#' @export
plotPlate <- function(dat, plotCol = "SampleGroup", wellCols = NULL){
    nCols <- ceiling(nrow(dat)/8)
    plt <- dat %>%  
        mutate_at(plotCol, as.character) %>%  
        ggplot(aes(x=Column, y=RowID))+
        geom_tile(aes_string(fill=plotCol), colour="black")  +
        scale_x_continuous(breaks=1:nCols) +
        scale_y_reverse(labels=LETTERS[1:8], breaks=c(1:8)) +
        theme(axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank()) +
        ggtitle(plotCol) +
        labs(fill=plotCol)
    if(is.null(wellCols)){ wellCols <- getCols(dat, plotCol) }
    plt <- plt + scale_fill_manual(values = wellCols)
    if(any(str_detect(colnames(dat), "[Rr]eplicate"))){
        plt <- plt + geom_text(aes(label=Replicate), size=7)
    }
    plt
}

outputPlatePlot <- function(batchColumn, datTab, outFileName){
    outFileName  <- str_c(outFileName, ".", batchColumn, ".png")
    p1 <- plotPlate(datTab, batchColumn)
    ggsave(outFileName, plot = p1)
} 

################################################################################


#' Fix metadata template
#'
#' Load and transform our CRUK CI Bioinfomratics Core metadata template into a
#' table suitable for the PlateLayout package
#' 
#' @param xlsFile The metadata template Excel file
#' @details The function takes the path for an Excel file containing completed 
#' CRUK CI metadata template and transforms it into a table suitable for use
#' with the \code{\link{randomizeSinglePlate}} function. 
#' @details When loading the Excel file the function skips the first row, 
#' making the headers from the second (hidden) row in the template. The third
#' row in the template contains a contents description for each column and so
#' is removed from the table.
#' @details Any columns that are completely empty are removed. \code{NA} is
#' replaced with a blank characer vector.
#' @return The cleaned metadata as a tibble and exports a tsv version
#' @examples
#' metadataFile <- system.file("extdata", "metadata_template_example.xls",
#'                            package = "PlateLayout")
#' metadataTable <- loadMetaForm(metadataFile)
#' metadataTable
#' @export
loadMetaForm <- function(xlsFile){
  outNam <- xlsFile %>% 
    str_replace_all(" ", "_") %>% 
    str_replace("\\.[xls]+$", ".tsv")
  message("Transforming '", xlsFile, "' to '", outNam, "'")
  dat <- read_excel(xlsFile, skip=1) %>%  
      slice(-1) %>%  
      mutate_all(as.character) %>%  
      rename(SampleGroup="Group") %>% 
      select_if(~any(!is.na(.x))) %>%  
      mutate_all(~ifelse(is.na(.x), "", .x))
  write_tsv(dat, outNam)
  dat
}

