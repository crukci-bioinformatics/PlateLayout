
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

# plot the plate with wells coloured according to a column in the samplesheet
# if the column is the "SampleGroup" the replicate numbers are also printed
library(RColorBrewer)

plotPlate <- function(dat, plotCol, boxColCodes = NULL){
    textSize <- ifelse(plotCol=="SampleGroup", 7, 0)
    nCols <- ceiling(nrow(dat)/8)
    len <- pull(dat, plotCol) %>% unique() %>% length()
    if(is.null(boxColCodes)) boxColCodes <- brewer.pal(len, "Set3")
    ggplot(dat, aes(x=Column, y=RowID))+
        geom_tile(aes_string(fill=plotCol), colour="black")  +
        geom_text(aes(label=Replicate), size=textSize) +
        scale_x_continuous(breaks=1:nCols) +
        scale_y_reverse(labels=LETTERS[1:8], breaks=c(1:8)) +
        scale_fill_manual(values=boxColCodes) +
        theme(axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank()) +
        ggtitle(plotCol) +
        labs(fill=plotCol)
}

savePlots <- function(batchColumn, datTab, outFileName){
    outFileName  <- str_c(outFileName, ".", batchColumn, ".png")
    p1 <- plotPlate(datTab, batchColumn)
    ggsave(outFileName, plot = p1)
} 
