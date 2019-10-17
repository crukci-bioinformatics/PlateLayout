randomizeSinglePlate <- function(designSheet, 
                                 outputFile, 
                                 primaryGroup = "SampleGroup", 
                                 batchColumns = NULL, 
                                 nIter = 10000,
                                 nCores = 4){
    # Set parallel processing
    plan(multiprocess, workers = nCores)

    # read sample sheet
    samsht <- read_tsv(designSheet, col_types = cols(.default = "c")) 

    # Modify batch columns
    bCols <- batchColumns %>%  
        str_split(",", simplify = TRUE) %>%  
        c(primaryGroup) %>%  
        str_subset("Replicate", negate = TRUE) %>%  
        unique()

    message("Run ", nIter, " random layouts")
    layouts <- future_map(1:nIter, 
                          ~plateRandomisation(samsht, primaryGroup = primaryGroup),
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
