# This script takes an experimental design table and generates a pseudo-block-randomised plate layout for it. It can cope with multiple plates - it will randomly spread the samples evenly between the plates. 
# The input table should be a csv file with a minimum of three columns headed "SampleName", SampleGroup" and "Replicate". An optional "Batch" column may also be specified. Any additional columns will be ignored.
# The outputs are one image file showing the plate layout, indicating both sample group and replicate number, and a table with the layout. The output table will include all the data in the input table, plus 3 or 4 additional columns indicating the plate number (if multiple plates), the row of the well (A-H), the column of the well (1-12) and the coordinate of the well (A1, A2...H7,H8)

library(optparse)

#Input options
options_list <- list(
  make_option(c('--layoutSheet', '-l'), type='character', help="Path to plate layout file - required", dest="layoutSheet", default="<undefined>"),
  make_option(c('--output', '-o'), type='character', help="Output Filename - optional", dest="outputFile", default="<undefined>"),
  make_option(c('--batchColumns', '-b'), type='character', help="Column to be plotted in multiplate experiments", dest="batCol", default="<undefined>")
)

#read options
parser <- OptionParser(option_list=options_list, usage="%prog [options]")
opts <- parse_args(parser)
layoutFile <- opts$layoutSheet
outputFile <- opts$outputFile
batchColumns <- opts$batCol

# set options that that have not beem provided
if(outputFile == "<undefined>") { outputFile <- gsub("[[:alnum:]]*$", "PlateLayoutPlots", basename(layoutFile)) }

######################################################################################################################################################################
# setwd("/run/user/1952417101/gvfs/sftp:host=clust1-headnode-1/home/sawle01/Groups/Ponder/20170109_OreillyM_PB_miRNAseq/PlateLayout/")
# layoutFile <- "miRNAseq.PlateLayout.csv"
# outputFile <- "TestRun1"
# batchColumns <- "Tissue,Outcome,Condition,SmokingStatus,Sex"


#########################################################################################################
## FUNCTIONS
#########################################################################################################

## Plate layout plot - shows the wells as a grid with sample groups represented by the colour and the replicate number printed in each well
testPlot <- function(plotDat, conditionToPlot, boxColCodes){
  textSize <- 0
  if(conditionToPlot=="SampleGroup") textSize <- 7
  ggplot(plotDat)+
    geom_tile(aes(x=Column, y=Row, fill=plotDat[[conditionToPlot]]), colour="black")  +
    geom_text(aes(x=plotDat$Column, y=plotDat$Row, label=ifelse(is.na(plotDat$Replicate), "", plotDat$Replicate)), size=textSize) +
    scale_x_continuous(breaks=1:NumberOfColumns) +
    scale_y_reverse(labels=LETTERS[1:8], breaks=c(1:8)) +
    scale_fill_manual(breaks=levels(plotDat[[conditionToPlot]]), values=boxColCodes) +
    theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
          axis.title.x=element_blank(), axis.title.y=element_blank()) +
    ggtitle(conditionToPlot) +
    labs(fill=conditionToPlot)
}

## Plot multi-plate bar plot for assessing distribution of factors between the plates
#plotDat <- tabout
multiplatePlot <- function(plotDat, batches){
  
  #remove the genomics controls and water, and transform the plate numbers in to a factor
  plotDat <- plotDat %>%  
    filter(!SampleName%in%c("Water", "GenomicsControl")) %>% 
    mutate(PlateNumber=factor(PlateNumber))
  
  #plot the distribution of sample groups across plates
  sampleGroupPlot <- ggplot(plotDat)+
    geom_bar(aes(x=SampleGroup, fill=PlateNumber), position="dodge")+
    theme(axis.text.x=element_text(angle=90))+
    ggtitle("SampleGroup")
  print(sampleGroupPlot)
  
  #plot the distribution of replicate groups across plates
  replicatePlot <- ggplot(plotDat)+
    geom_bar(aes(x=Replicate, fill=PlateNumber), position="dodge")+
    theme(axis.text.x=element_text(angle=90))+
    ggtitle("Replicate")
  print(replicatePlot)
  
  #plot the distribution of other requested factors across plates
  if(batches!="<undefined>"){
    allbatches <- strsplit(batches, ",")[[1]]
    for(colName in allbatches){
      colNamePlot <- ggplot(plotDat, aes_string(x=colName, fill="PlateNumber"))+
        geom_bar(position="dodge")+
        theme(axis.text.x=element_text(angle=90))+
        ggtitle(colName)
      print(colNamePlot)
    }
  }
}
#########################################################################################################
## MAIN SCRIPT
#########################################################################################################

library(dplyr)
library(ggplot2)
library(RColorBrewer)


dat <- read.csv(layoutFile)
NumberOfPlates <- max(c(1, dat$PlateNumber))
groupsList <- levels(dat$SampleGroup)
groupsList <- groupsList[order(groupsList)]

#fix up the factors in the batch columns
if(batchColumns!="<undefined>"){
    allbatches <- strsplit(batchColumns, ",")[[1]]
    for(colName in allbatches){
        newLevels <- c(levels(dat[,colName])[!levels(dat[,colName])%in%c("Water", "GenomicsControl")], "Water", "GenomicsControl")
        dat[,colName] <- factor(dat[,colName], levels = newLevels)
    }
}
  
# randomise each plate
plateID <- ""
for(thisPlateNumber in 1:NumberOfPlates){
    
    if(NumberOfPlates>1){
        plateID <- paste(".Plate_", thisPlateNumber, sep="")
        #subset the sample data for this plate
        plateDat <- filter(dat, PlateNumber==thisPlateNumber)
    }else{
        plateDat <- dat
    }

    #replace the well row letters with numbers
    plateDat <- plateDat %>% 
        mutate(Row=match(Row, LETTERS))
    
    NumberOfColumns <- max(plateDat$Column)
    nWater <- length(which(plateDat$SampleName=="Water"))
    nGenomicControls <- length(which(plateDat$SampleName=="GenomicsControl"))
    
    # plot the plate layout and save
     #generate a colour coding for plotting to keep consistency across different plates
    if(thisPlateNumber==1){
        colCodes <- colorRampPalette(brewer.pal(12,"Set3"))(length(levels(plateDat$SampleGroup))-2)
        if(nWater>0) colCodes <- c(colCodes, "#0000FF")
        if(nGenomicControls>0) colCodes <- c(colCodes, "#000000")
    }
    #plot plate layout
    testPlot(plateDat, "SampleGroup", colCodes)
    plotWidth <- (1.1*NumberOfColumns)+(max(nchar(groupsList))*0.075)+0.15
    ggsave(paste(outputFile, plateID, ".png", sep=""), height = 5, width=plotWidth)  
    
    if(batchColumns!="<undefined>"){
      allbatches <- strsplit(batchColumns, ",")[[1]]
      #generate a colour coding for plotting each batch to keep consistency across different plates
      if(thisPlateNumber==1){ 
        colCodesBatches <- vector("list", length(allbatches)) 
        names(colCodesBatches) <- allbatches
        }
      for(colName in allbatches){
        if(thisPlateNumber==1){
          tempCols <- colorRampPalette(brewer.pal(12,"Set3"))(length(levels(plateDat[[colName]]))-2)
          if(nWater>0) tempCols <- c(tempCols, "#0000FF")
          if(nGenomicControls>0) tempCols <- c(tempCols, "#000000")
          colCodesBatches[[colName]] <- tempCols
        }
        #plot
        testPlot(plateDat, colName, colCodesBatches[[colName]])
        plotWidth <- (1.1*NumberOfColumns)+1.5
        ggsave(paste(outputFile, plateID, ".", colName, ".png", sep=""), height = 5, width=plotWidth)  
      }
    }
}

#create a pdf with plots showing the distribution of the different factors across plates
if(NumberOfPlates>1){
    pdf(paste(outputFile, ".Multiplate_Distribution_Check.pdf", sep=""))
    multiplatePlot(dat[dat$SampleName!="GenomicsControl",], batchColumns)
    dev.off()
}