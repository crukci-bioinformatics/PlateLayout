# This script takes an experimental design table and generates a pseudo-block-randomised plate layout for it. It can cope with multiple plates - it will randomly spread the samples evenly between the plates. 
# The input table should be a csv file with a minimum of three columns headed "SampleName", SampleGroup" and "Replicate". An optional "Batch" column may also be specified. Any additional columns will be ignored.
# The outputs are one image file showing the plate layout, indicating both sample group and replicate number, and a table with the layout. The output table will include all the data in the input table, plus 3 or 4 additional columns indicating the plate number (if multiple plates), the row of the well (A-H), the column of the well (1-12) and the coordinate of the well (A1, A2...H7,H8)

library(optparse)

#Input options
options_list <- list(
  make_option(c('--designSheet', '-d'), type='character', help="Path to experimental design file - required", dest="designSheet", default="<undefined>"),
  make_option(c('--output', '-o'), type='character', help="Output Filename - optional", dest="outputFile", default="<undefined>"),
  make_option(c('--batchColumns', '-b'), type='character', help="Column to be plotted in multiplate experiments", dest="batCol", default="<undefined>"),
  make_option(c('--runNumber', '-r'), type='integer', help="Number of times to run the script", dest="runNum", default="1"),
  make_option(c('--noGenomicsControls', '-G'), action="store_false", default=TRUE, help="Do not add Genomics controls (generally only used when the lab is doing the library prep rather than Genomics)", dest="genCtrls")
)

#read options
parser <- OptionParser(option_list=options_list, usage="%prog [options]")
opts <- parse_args(parser)
designFile <- opts$designSheet
outputFile <- opts$outputFile
batchColumns <- opts$batCol
runNumber <- opts$runNum
genCtrls <- opts$genCtrls

# set options that that have not beem provided
if(outputFile == "<undefined>") { outputFile <- gsub("[[:alnum:]]*$", "PlateLayout", basename(designFile)) }

#########################################################################################################
## FUNCTIONS
#########################################################################################################

# Randomise order of rows or columns
randBlock <- function(x){
  randomOrder=sample(1:max(x))
  randomOrder[x]
}

## Plate layout plot - shows the wells as a grid with sample groups represented by the colour and the replicate number printed in each well
testPlot <- function(plotDat, conditionToPlot, boxColCodes){
  textSize <- 0
  if(conditionToPlot=="SampleGroup") textSize <- 7
  ggplot(plotDat)+
    geom_tile(aes(x=Column, y=Row, fill=plotDat[[conditionToPlot]]), colour="black")  +
    geom_text(aes(x=plotDat$Column, y=plotDat$Row, label=plotDat$Replicate), size=textSize) +
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
    filter(!SampleGroup%in%c("Water", "GenomicsControl")) %>% 
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
options(warn=2)

library(dplyr)
library(ggplot2)
library(RColorBrewer)


runID <- ""
# run the script as many times as requested
for(attemptNumber in 1:runNumber){
  if(runNumber>1) {
    cat(paste("Run Number:", attemptNumber, "\n"))
    runID <- paste(".Run", attemptNumber, sep="") 
  }
  
  # Get sample data
  dat <- read.csv(designFile, stringsAsFactors=F) %>% 
    tbl_df()
  
  #Set the number Wells on a plate
  if(genCtrls){
    WellsInOnePlate <- 94 # 2 wells are reserved for the Genomics controls
  }else{
    WellsInOnePlate <- 96
  }
  
  # Calculate number of plates required and number of samples per plate. The aim here is to distribute samples as evenly as possible between multiple plates, but to uses the minimum number of columns. i.e. e.g if there are 100 samples, splitting the samples evenly would put 50 samples on each plate. This would use 7 columns on each plate with 6 empty wells on each plate. More efficiently we could put 48 samples on one plate, using 6 full columns, and 52 samples on the other using 7 columns with 4 empty wells.
  NumberOfSamples <- nrow(dat)
  noS <- NumberOfSamples
  # the following loop optimises the number columns used across the plates and outputs a vector of the number of columns on each plate. If you can figure out a more elegant way of doing this please let me know.
  count12s <- 0
  count12sold <- -1
  while(count12sold!=count12s){
    noS <- noS+(count12s*2)
    NumberOfColumns <- ceiling(noS/8)
    NumberOfPlates <- ceiling(NumberOfColumns/12)
    columnsPerPlate <- NumberOfColumns%/%NumberOfPlates
    extraColumns <- NumberOfColumns%%NumberOfPlates
    columnsOnEachPlate <- rep(columnsPerPlate, NumberOfPlates)
    if(extraColumns>0) columnsOnEachPlate[1:extraColumns] <- columnsPerPlate+1
    count12sold <- count12s
    if(genCtrls) count12s <- length(which(columnsOnEachPlate==12)) #this will cause the will cause the loop to cycle, adding two samples for each complete 12 column plate to represent the genomics controls 
  }
  #create a vector of the number of samples that can go on each plate given the number of columns - the max is 94 if genomics controls are to be added
  samplesOnEachPlate <- columnsOnEachPlate*8
  samplesOnEachPlate[samplesOnEachPlate>WellsInOnePlate] <- WellsInOnePlate
  #there may be less samples that the total number of available wells. Remove the unecessary wells from plates with the larger number of columns
  nLess <- sum(samplesOnEachPlate)-nrow(dat)
  if(nLess>0){
    whi <- which(samplesOnEachPlate==max(samplesOnEachPlate))
    samplesOnEachPlate[whi] <- samplesOnEachPlate[whi]-(nLess%/%length(whi))
    extras <- sample(whi, nLess%%length(whi))
    samplesOnEachPlate[extras] <- samplesOnEachPlate[extras]-1
  }
  #generate a vector to bind to the sample table that will assign each sample to a plate. This is organised to distribute the sample groups as evenly as possible across the plates
  samplePlateAssignment <- vector()
  minNSamples <- min(samplesOnEachPlate)
  for(i in 1:length(samplesOnEachPlate)){
    if(samplesOnEachPlate[i]==minNSamples){
      samOrd <- 1:minNSamples
    }else{
      samOrd <- c(1:minNSamples , sample(1:minNSamples, 8))[1:samplesOnEachPlate[i]] # otherwise we'd end up with just the larger plates at the end of the list meaning that some sample groups could potentially end up on just one plate
    }
    temp <- cbind(rep(i, samplesOnEachPlate[i]), samOrd)
    samplePlateAssignment <- rbind(samplePlateAssignment, temp)
  }
  samplePlateAssignment <- samplePlateAssignment[order(samplePlateAssignment[,2], samplePlateAssignment[,1]),]
  
  
  ###################
  #order the factor level of the sample groups. Add water and genomics controls to the factor levels for compatibility with the data frames for these that we may bind later. Also change the replicate column to character.
  groupsList <- unique(dat$SampleGroup)
  groupsList <- groupsList[order(groupsList)]
  dat <- mutate(dat, SampleGroup=factor(SampleGroup, levels=c(groupsList, "Water", "GenomicsControl"))) %>% 
      mutate(Replicate=as.character(Replicate))
      #mutate(Replicate=factor(Replicate, levels=c(unique(plateDat$Replicate), "")))
  
  #change any additional batch columns into factors including levels for the water and genomics controls which may be added later
  if(batchColumns!="<undefined>"){
    allbatches <- strsplit(batchColumns, ",")[[1]]
    for(colName in allbatches){
      dat[,colName] <- as.factor(unlist(dat[,colName]))
      levels(dat[[colName]]) <- c(levels(dat[[colName]]),  "Water", "GenomicsControl")
    }
  }
  
  # Randomise the order of the Groups and the samples and then assign to plates
  Groups <- unique(dat$SampleGroup)
  OrderOfGroups <- sample(1:length(Groups))
  dat <- dat %>% 
    mutate(GroupNumber=OrderOfGroups[match(SampleGroup, Groups)]) %>% 
    mutate(SampleNumber=sample(1:NumberOfSamples)) %>% 
    arrange(GroupNumber, SampleNumber) %>% 
    mutate(PlateNumber=samplePlateAssignment[1:nrow(dat),1])
  
  # create a data frame to contain all plates for final the output
  tabout <- tbl_df(data.frame())
  
  # randomise each plate
  plateID <- ""
  for(thisPlateNumber in 1:NumberOfPlates){
    
    if(NumberOfPlates>1){ plateID <- paste(".Plate_", thisPlateNumber, sep="") }
    
    #subset the sample data for this plate
    plateDat <- filter(dat, PlateNumber==thisPlateNumber)
    
    # calculate how many water controls and genomics controls to add
    # We must fill up complete columns
    # up to 2 empty wells will be filled with Genomics controls, any remaining wells are water
    nSamples <- nrow(plateDat)
    nEmptyWells <- 8-(nrow(plateDat)%%8)
    if(nEmptyWells==8) nEmptyWells <- 0
    NumberOfColumns <- ceiling(nSamples/8)
    if(genCtrls){
      nGenomicControls <- min(c(2, nEmptyWells)) # Fill any empty wells with up to 2 Genomic Controls
      nWater <- nEmptyWells - nGenomicControls
    }else{
      nGenomicControls <- 0
      nWater <- nEmptyWells
    }
    
    # create a matrix for the wells that will be used, moving diagonally across the plate 
    AllWells <- data.frame(Column=0, Row=rep(1:8, NumberOfColumns))
    COL <- 0
    for(j in 1:nrow(AllWells)){
      COL <- COL + 1
      if(paste(COL, AllWells$Row[j])%in%paste(AllWells$Column[1:(j-1)], AllWells$Row[1:(j-1)])) { COL <- COL + 1 }
      AllWells$Column[j] <- COL
      if(COL==NumberOfColumns) COL <- 0
    }
    
    #Add the WaterSamples. Water added into random positions within the data frame
    if(nWater>0){
      # create random indexes for the positions for the water wells to be inserted into the plate table and add 0.75. When the two dataframe are bound and sorted by "SampleIndex" this will cause the water wells to be inserted into the rows after the random rows genereate
      # i.e. if sample gives us rows 5 and 13, the water indexes will be 5.75 and 13.75, meaning that after sorting they will be in rows 6 and 14
      waterWells <- rep(sample(0:nSamples, 1)+0.75, nWater)
      
      #Generate a dataframe for the water controls
      pseudoSamples <- rep("Water", nWater)
      pseudoSamples <- factor(pseudoSamples, levels=c(groupsList, "Water", "GenomicsControl"))
      pseudoSamplesTab <- tbl_df(data.frame(SampleGroup=pseudoSamples, 
                                         SampleName=pseudoSamples, 
                                         SampleIndex=waterWells,
                                         Replicate=as.character(rep("", nWater)),
                                         #Replicate=factor(rep("", nWater), levels=c(unique(plateDat$Replicate), "")),
                                         PlateNumber=rep(thisPlateNumber, nWater)
                              ))
      #add batch columns if necessary
      if(batchColumns!="<undefined>"){
        allbatches <- strsplit(batchColumns, ",")[[1]]
        temp <- as.data.frame(matrix("Water", nr=nrow(pseudoSamplesTab), nc=length(allbatches)))
        for(colName in 1:ncol(temp)){temp[,colName] <- factor(temp[,colName], levels=levels(plateDat[[match(allbatches[colName], colnames(plateDat))]]))}
        colnames(temp) <- allbatches
        pseudoSamplesTab <- cbind(pseudoSamplesTab, temp)
      }
      #bind and sort the two data frames
      plateDat <- plateDat %>% 
        mutate(SampleIndex=1:nSamples) %>% 
        bind_rows(pseudoSamplesTab) %>% 
        arrange(SampleIndex) %>% 
        select(-SampleIndex) 
    }
    
    #Add the Genomic Controls. GC are added in the rows corresponding to the final two wells used on the plate (they will move during randomisation but when we swap them back this ensures the correct relationship for the two samples they are swapped with - i.e. they were in the same column - the only thing we can know about them at the moment)
    if(nGenomicControls>0){

      # find the row index in the array for the bottom wells in the last column for the genomics controls to be placed into
      GCWells <- c(which(AllWells$Column==NumberOfColumns&AllWells$Row==8), which(AllWells$Column==NumberOfColumns&AllWells$Row==7))[1:nGenomicControls]-0.5
      if(nGenomicControls==2) GCWells[which.max(GCWells)] <- GCWells[which.max(GCWells)]-1
      
      #Generate a dataframe for the Genomics Controls
      pseudoSamples <-rep("GenomicsControl", nGenomicControls)
      pseudoSamples <- factor(pseudoSamples, levels=c(groupsList, "Water", "GenomicsControl"))
      pseudoSamplesTab <- tbl_df(data.frame(SampleGroup=pseudoSamples, 
                                         SampleName=pseudoSamples, 
                                         SampleIndex=GCWells,
                                         Replicate=rep("", nGenomicControls),
                                         #Replicate=factor(rep("", nGenomicControls), levels=c(unique(plateDat$Replicate, ""))),
                                         PlateNumber=rep(thisPlateNumber, nGenomicControls)
      ))
      #add batch columns if necessary
      if(batchColumns!="<undefined>"){
        allbatches <- strsplit(batchColumns, ",")[[1]]
        temp <- as.data.frame(matrix("GenomicsControl", nr=nrow(pseudoSamplesTab), nc=length(allbatches)))
        for(colName in 1:ncol(temp)){temp[,colName] <- factor(temp[,colName], levels=levels(plateDat[[match(allbatches[colName], colnames(plateDat))]]))}
        colnames(temp) <- allbatches
        pseudoSamplesTab <- cbind(pseudoSamplesTab, temp)
      }
      #bind and sort the two data frames
      plateDat <- plateDat %>% 
        mutate(SampleIndex=1:(nSamples+nWater)) %>% 
        bind_rows(pseudoSamplesTab) %>% 
        arrange(SampleIndex) %>% 
        select(-SampleIndex) 
    }
    
    
    #Now that we have a complete sample table we can add the well coordinates - this is done systematically to create a pseudo-blocking effect, we will randomise afterwards
    plateGroupList <- unique(c(groupsList, rep("Water", nWater), rep("GenomicsControls", nGenomicControls)))
    plateDat <- plateDat %>% 
      bind_cols(AllWells) #%>% 
      #droplevels
    
    #Randomise rows and then columns
    plateDat <- plateDat %>% 
      mutate(Row=randBlock(Row)) %>% 
      mutate(Column=randBlock(Column))
    
    #finally switch the Genomics controls with the end wells if necessary
    if(nGenomicControls>0){
      GCRows <- which(plateDat$Column==NumberOfColumns&plateDat$Row>(8-nGenomicControls))
      AreGC <- which(plateDat$SampleName=="GenomicsControl")
      #if any are already in the right place, remove them
      GCrows <- GCRows[!GCRows%in%AreGC]
      areGC <- AreGC[!AreGC%in%GCRows]
      #adjust the row and columns accordingly
      plateDat <- plateDat %>% 
        mutate(NewRow=Row) %>% 
        mutate(NewColumn=Column) %>% 
        mutate(NewRow=replace(NewRow, GCrows, plateDat$Row[areGC]))    %>% 
        mutate(NewRow=replace(NewRow, areGC, plateDat$Row[GCrows])) %>% 
        mutate(NewColumn=replace(NewColumn, GCrows, plateDat$Column[areGC]))    %>% 
        mutate(NewColumn=replace(NewColumn, areGC, plateDat$Column[GCrows])) %>% 
        select(-Row, -Column) %>% 
        rename(Row=NewRow) %>% 
        rename(Column=NewColumn)
    }

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
    ggsave(paste(outputFile, runID, plateID, ".png", sep=""), height = 5, width=plotWidth)  
    
    if(batchColumns!="<undefined>"){
      allbatches <- strsplit(batchColumns, ",")[[1]]
      #generate a colour coding for plotting each batch to keep consistency across different plates
      if(thisPlateNumber==1){ 
        colCodesBatches <- vector("list", length(allbatches)) 
        names(colCodesBatches) <- allbatches
        }
      for(colName in allbatches){
        if(thisPlateNumber==1){
          colCodesBatches[[colName]] <- colorRampPalette(brewer.pal(12,"Set3"))(length(levels(plateDat[[colName]]))-2)
          if(nWater>0) colCodesBatches[[colName]] <- c(colCodesBatches[[colName]], "#0000FF")
          if(nGenomicControls>0) colCodesBatches[[colName]] <- c(colCodesBatches[[colName]], "#000000")
        }
        #plot
        testPlot(plateDat, colName, colCodesBatches[[colName]])
        plotWidth <- (1.1*NumberOfColumns)+1.5
        ggsave(paste(outputFile, runID, plateID, ".", colName, ".png", sep=""), height = 5, width=plotWidth)  
      }
    }
    
    #replace the welll row numbers with letters and create the well index column
    plateDat <- plateDat %>% 
      mutate(Row=LETTERS[Row]) %>% 
      mutate(Well=paste(Row, Column, sep=""))
    
    # bind the plate to the output table
    tabout <- plateDat %>% 
      select(-GroupNumber, -SampleNumber) %>% 
      mutate(Well=paste(Row, Column, sep="")) %>% 
      bind_rows(tabout)
  }
  
  # if there is only one plate remove the plate number column. If there are multiple plates create a pdf with the bar plots for the distribution of factors between the plates.
  if(NumberOfPlates==1){ 
    dat <- select(dat, -PlateNumber)
  }else{
    pdf(paste(outputFile, ".Multiplate_Distribution_Check", runID, ".pdf", sep=""))
     multiplatePlot(tabout, batchColumns)
    dev.off()
  }
  tabout <- tabout %>% 
    arrange(PlateNumber, Column, Row)
  #write the final output table
  write.csv(tabout, paste(outputFile, runID, ".csv", sep=""), row.names = F, quote=F)
}
