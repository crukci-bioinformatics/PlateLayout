# This script takes an experimental design table and generates a pseudo-block-randomised plate layout for it. It can cope with multiple plates - it will randomly spread the samples evenly between the plates. 
# The input table should be a csv file with a minimum of three columns headed "SampleName", SampleGroup" and "Replicate". Any additional columns will be ignored.
# If the samples are to be spread across multiple plates you can include a "PlateNumber" column, and manually assign the samples to specific plates.
# The outputs are one image file showing the plate layout, indicating both sample group and replicate number, and a table with the layout. The output table will include all the data in the input table, plus 3 or 4 additional columns indicating the plate number (if multiple plates), the row of the well (A-H), the column of the well (1-12) and the coordinate of the well (A1, A2...H7,H8)

library(optparse)

#Input options
options_list <- list(
    make_option(c('--designSheet', '-d'), type='character', help="Path to experimental design file - required", dest="designSheet", default="<undefined>"),
    make_option(c('--output', '-o'), type='character', help="Output Filename - optional", dest="outputFile", default="<undefined>"),
    make_option(c('--batchColumns', '-b'), type='character', help="Column to be plotted in multiplate experiments", dest="batCol", default="<undefined>"),
    make_option(c('--maxWells', '-w'), type='integer', help="The maximum number of wells to use per plate", dest="maxWells", default="96"),
    make_option(c('--noGenomicsControls', '-G'), action="store_false", default=TRUE, help="Do not add Genomics controls (generally only used when the lab is doing the library prep rather than Genomics)", dest="genCtrls")
)

#read options
parser <- OptionParser(option_list=options_list, usage="%prog [options]")
opts <- parse_args(parser)
designFile <- opts$designSheet
outputFile <- opts$outputFile
BatchColumns <- opts$batCol
WellsInOnePlate <- opts$maxWells
genCtrls <- opts$genCtrls

# set options that that have not beem provided
if(outputFile == "<undefined>") { outputFile <- gsub("[[:alnum:]]*$", "PlateLayout", basename(designFile)) }

#########################################################################################################
## FUNCTIONS
#########################################################################################################

# determine the number of plates to use and the number of columns on each plate. see NOTE 1 at end of script
getNumberOfPlates <- function(numberOfSamples, columnsInOnePlate){
    noS <- numberOfSamples
    count12s <- 0
    count12sold <- -1
    # this loop determines the number of columns to use one each plate, accounting for two genomic controls in full plates. If you can determine a more elegant way of doing this please let me know.
    while(count12sold!=count12s){
        noS <- noS+(count12s*2)
        numberOfColumns <- ceiling(noS/8)
        numberOfPlates <- ceiling(numberOfColumns/columnsInOnePlate)
        columnsPerPlate <- numberOfColumns%/%numberOfPlates
        extraColumns <- numberOfColumns%%numberOfPlates
        columnsOnEachPlate <- rep(columnsPerPlate, numberOfPlates)
        columnsOnEachPlate[seq(0, extraColumns)] <- columnsOnEachPlate[seq(0, extraColumns)]+1
        count12sold <- count12s
        if(genCtrls) count12s <- length(which(columnsOnEachPlate==12)) #add two samples for each complete 12 column plate to represent the genomics controls and recalculate if necessary
    }
    return(columnsOnEachPlate)
}

# determine the number of samples to put on each plate
getSamplesPerPlate <- function(numberOfSamples, columnsOnEachPlate, wellsInOnePlate){
    numberOfPlates <- length(columnsOnEachPlate)
    #create a vector of the number of available wells on each plate given the number of columns - the max is 94 if genomics controls are to be added
    samplesOnEachPlate <- sapply(columnsOnEachPlate, function(x){min(x*8, wellsInOnePlate)})
    
    #there may be less samples that the total number of available wells.
    nLess <- sum(samplesOnEachPlate)-numberOfSamples
    if(nLess>0){
        remWell <- rep(nLess%/%numberOfPlates, numberOfPlates)
        extras <- nLess%%numberOfPlates
        remWell[seq(0, extras)] <- remWell[seq(0, extras)] + 1
        samplesOnEachPlate <- samplesOnEachPlate-remWell
    }
    return(samplesOnEachPlate)
}

# determine the distribution of samples across multiple plates
distributeSamples <- function(datTab, batchColumns, samplesOnEachPlate){

    # The algorithm to assign the samples to multiple plates contains two random elements. 
    # We will run it multiple times and assess which results in the most even distribution of the "batch" factors
    repTests <- 10000 #The number of test we will do
    # names of columns to test
    testColumns <- "SampleGroup" # we will always assess at least the SampleGroup column
    if(batchColumns!="<undefined>"){
        testColumns <- c(testColumns, strsplit(batchColumns, ",")[[1]])
    }
    # output matrix of the distribution statitics for each factor for each test
    distrFactors <- matrix(nc=length(testColumns), nr=repTests)
    colnames(distrFactors) <- testColumns
    # lists to contain the random orders generated
    listOrderOfGroups <- list()
    listOrderOfSamples <- list()
    listSamplePlateAssignment <- list()
    for(i in 1:repTests){
        # Generate a vector to bind to the sample table that will assign each sample to a plate. 
        # This is organised to distribute the sample groups as evenly as possible across the plates.
        minSam <- min(samplesOnEachPlate)
        samplePlateAssignment <- rep(1:length(samplesOnEachPlate), minSam)
        samDiff <- samplesOnEachPlate-minSam
        addVec <- as.vector(sapply(which(samDiff>0), function(x){rep(x, samDiff[x])}))
        samplePlateAssignment <- c(samplePlateAssignment, addVec)[order(c(1:length(samplePlateAssignment), sample(0:length(samplePlateAssignment), length(addVec))+0.5))]
        
        # Randomise the order of the Groups and the samples and then assign to plates
        Groups <- unique(datTab$SampleGroup)
        orderOfGroups <- sample(1:length(Groups))
        orderOfSamples <- sample(1:NumberOfSamples)
        tempDat <- datTab %>% 
            mutate(GroupNumber=orderOfGroups[match(SampleGroup, Groups)]) %>% 
            mutate(SampleNumber=orderOfSamples) %>% 
            arrange(GroupNumber, SampleNumber) %>% 
            mutate(PlateNumber=samplePlateAssignment)
        
        # See how far the distrtibution of each factor differs from ideal
        for(thisBatch in testColumns){
            distrFactors[i,thisBatch] <- getDistrFactor(tempDat, thisBatch)
        }
        listOrderOfGroups[[i]] <- orderOfGroups
        listOrderOfSamples[[i]] <- orderOfSamples
        listSamplePlateAssignment[[i]] <- samplePlateAssignment
    }
    
    # Select the optimum distribution
    minFactors <- apply(distrFactors, 2, min)
    factorDiff <- apply(distrFactors, 1, function(x){sum(abs(x-minFactors))})
    write.table(distrFactors, "PlateDistributionFactors.tab", sep="\t", col.names=T, row.names=F, quote=F)
    optDistr <- which(factorDiff==min(factorDiff))[1]
    # assign plate numbers
    datTab <- datTab %>% 
        mutate(GroupNumber=listOrderOfGroups[[optDistr]][match(SampleGroup, Groups)]) %>% 
        mutate(SampleNumber=listOrderOfSamples[[optDistr]]) %>% 
        arrange(GroupNumber, SampleNumber) %>% 
        mutate(PlateNumber=listSamplePlateAssignment[[optDistr]])
    return(datTab)
}

# test the distribution of "batch" factors across plates
# The returned statistic is the number of samples that would have to be moved to achieve the "optimum" distribution
# for all groups in the factor
getDistrFactor <- function(datTab, batchName){
    # list of different groups in the factor
    grpList <- unique(datTab[,batchName])
    nPlates <- max(datTab$PlateNumber)
    distrFactor <- 0
    for(grp in grpList){
        whi <- which(datTab[,batchName]==grp)
        nReps <- length(whi)
        # determine the optimum distribution of samples
        optDistr <- rep(nReps%/%nPlates, nPlates)
        optRem <- nReps%%nPlates
        if(optRem>0) optDistr[1:optRem] <- optDistr[1:optRem]+1
        # get the actual distribution
        actDistr <- aggregate(rep(1, nReps), by=list(datTab$PlateNumber[whi]), sum)[,2]
        actDistr <- c(actDistr, rep(0, length(optDistr)-length(actDistr)))
        actDistr <- sort(actDistr, decreasing = T)
        # determin the difference between the actual and optimum for this group
        distrDiff <- sum(abs(optDistr-actDistr))/2
        # add to the total statistic for the factor
        distrFactor <- distrFactor+distrDiff
    }
    return(distrFactor)
}

# Randomise order of rows or columns
randBlock <- function(x){
    randomOrder=sample(1:max(x))
    randomOrder[x]
}

# test the distribution of "batch" factors within a plate
# the base statistic is the count for each sample of the number of samples adjacent to it or something like that
getDistrScore <- function(datTab, batchName){
    datTab <- as.data.frame(datTab)
    distrScore <- 0
    for(thisRow in 1:nrow(datTab)){
        whi <- which(datTab[,batchName]==datTab[thisRow,batchName])
        myD <- abs(datTab$Column[whi]-datTab$Column[thisRow])+abs(datTab$Row[whi]-datTab$Row[thisRow])
        distrScore <- distrScore+length(which(myD<=2))+length(which(myD<=1))-2
    }
    return(distrScore)
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

#########################################################################################################
## MAIN SCRIPT
#########################################################################################################
options(warn=2)

library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Get sample data
dat <- read.csv(designFile, stringsAsFactors=F) %>% 
    tbl_df()

##Figure out the number of samples on each plate

#Set the number Wells & number of columns on a plate
ColumnsInOnePlate <- WellsInOnePlate/8
if(WellsInOnePlate%%8>0){
    stop("The script currently only works if the number of wells per plate is divisible by 8 i.e. using complete columns.\n You have specified ", WellsInOnePlate, " wells, which gives ", ColumnsInOnePlate)
}
if(WellsInOnePlate==96&genCtrls) WellsInOnePlate <- 94 # 2 wells are reserved for the Genomics controls if the plate is full

# Calculate number of plates required and number of samples per plate. 
#if a plate number column is included, check it, and if it is okay use this to calculate get these numbers
if("PlateNumber"%in%colnames(dat)){
    plateNumbers <- dat$PlateNumber
    SamplesOnEachPlate <- aggregate(plateNumbers, by=list(plateNumbers), length)[,2]
    if(any(SamplesOnEachPlate>WellsInOnePlate)) stop("One of the plates in the \"PlateNumber\" column contains more than the allowable number of wells (", WellsInOnePlate, ")")
    if(any(plateNumbers<1)) stop("One of the plates in the \"PlateNumber\" column is less than 1.")
    if(any(!is.integer(plateNumbers))) stop("The  \"PlateNumber\" column should only contain integers.")
    NumberOfSamples <- length(plateNumbers)
    ColumnsOnEachPlate <- ceiling(SamplesOnEachPlate/8)
    NumberOfPlates <- max(plateNumbers)
}else{
    NumberOfSamples <- nrow(dat)
    ColumnsOnEachPlate <- getNumberOfPlates(NumberOfSamples, ColumnsInOnePlate)
    NumberOfPlates <- length(ColumnsOnEachPlate)
    SamplesOnEachPlate <- getSamplesPerPlate(NumberOfSamples, ColumnsOnEachPlate, WellsInOnePlate)
}


# Get Number of samples on each plate

# Output results of above
message("Number of columns per plate: ", ColumnsInOnePlate)
message("Max number of wells per plate: ", WellsInOnePlate)
message("Total number of samples: ", NumberOfSamples)
message("Number of plates: ", NumberOfPlates)
message("Number of columns used on each plate: ", paste(ColumnsOnEachPlate, collapse=","))
message("Number of samples on each plate: ", paste(SamplesOnEachPlate, collapse=","))

#if the samples will fit on a single plate then assign all samples to "1", otherwise we need to optimise the distribution of the samples across the plates according to the different "batch" factors
if(!"PlateNumber"%in%colnames(dat)&NumberOfPlates==1){
    dat$PlateNumber <- 1
}else if(!"PlateNumber"%in%colnames(dat)&NumberOfPlates>1){
    #### distribute across multiple plates
    dat <- distributeSamples(dat, BatchColumns, SamplesOnEachPlate)
}

###################
#order the factor level of the sample groups. Add water and genomics controls to the factor levels for compatibility with the data frames for these that we may bind later. Also change the replicate column to character.
groupsList <- unique(dat$SampleGroup)
groupsList <- groupsList[order(groupsList)]
dat <- dat %>%
    mutate(SampleGroup=factor(SampleGroup, levels=c(groupsList, "Water", "GenomicsControl"))) %>%
    mutate(Replicate=as.character(Replicate))

#change any additional batch columns into factors including levels for the water and genomics controls which may be added later
if(BatchColumns!="<undefined>"){
    allbatches <- strsplit(BatchColumns, ",")[[1]]
    for(colName in allbatches){
        dat[,colName] <- as.factor(unlist(dat[,colName]))
        levels(dat[[colName]]) <- c(levels(dat[[colName]]),  "Water", "GenomicsControl")
    }
}
###################
message(dat$PlateNumber)

# create a data frame to contain all plates for final the output
tabout <- tbl_df(data.frame())

# randomise each plate
runNumber <- 10000
plateID <- ""
for(thisPlateNumber in 1:NumberOfPlates){
    
    if(NumberOfPlates>1) plateID <- paste(".Plate", thisPlateNumber, sep="_")
    # The plate layout algorithm first optimises for sample group and then randomises the rows and columns. In 
    # order to ensure the other "batch" factors are not spatially biased we will run this process a number of times
    # and then keep the one with the best distribution of samples
    testColumns <- "SampleGroup" # we will always assess at least the SampleGroup column
    if(BatchColumns!="<undefined>"){
        testColumns <- c(testColumns, strsplit(BatchColumns, ",")[[1]])
    }
    layoutScoreTab <- matrix(nc=length(testColumns), nr=runNumber)
    colnames(layoutScoreTab) <- testColumns
    layoutList <- list()
    
    for(thisRun in 1:runNumber){
        
        cat(paste("Plate", thisPlateNumber, "Run Number:", thisRun, "\n"))
        #subset the sample data for this plate
        plateDat <- filter(dat, PlateNumber==thisPlateNumber)
        
        #randomise the order of the groups and the replicates within each group
        Groups <- plateDat$SampleGroup
        orderOfGroups <- sample(1:length(Groups))
        orderOfSamples <- sample(1:nrow(plateDat))
        plateDat <- plateDat %>% 
            mutate(GroupNumber=orderOfGroups[match(SampleGroup, Groups)]) %>% 
            mutate(SampleNumber=orderOfSamples) %>% 
            arrange(GroupNumber, SampleNumber)
        
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
            pseudoSamplesTab <- tbl_df(data.frame(SampleGroup=factor("Water", levels=c(groupsList, "Water", "GenomicsControl")), 
                                                  SampleName="Water", 
                                                  SampleIndex=waterWells,
                                                  Replicate="",
                                                  PlateNumber=rep(thisPlateNumber, nWater), 
                                                  stringsAsFactors = F
            ))
            #add batch columns if necessary
            if(BatchColumns!="<undefined>"){
                allbatches <- strsplit(BatchColumns, ",")[[1]]
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
            pseudoSamplesTab <- tbl_df(data.frame(SampleGroup=factor("GenomicsControl", levels=c(groupsList, "Water", "GenomicsControl")),
                                                  SampleName="GenomicsControl", 
                                                  SampleIndex=GCWells,
                                                  Replicate="",
                                                  PlateNumber=rep(thisPlateNumber, nGenomicControls),
                                                  stringsAsFactors = F
            ))
            #add batch columns if necessary
            if(BatchColumns!="<undefined>"){
                allbatches <- strsplit(BatchColumns, ",")[[1]]
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
        
        # now that the layout is finalised we will test the distribution of the different factors and assign a score
        for(thisCol in testColumns){
            layoutScoreTab[thisRun,thisCol] <- getDistrScore(plateDat, thisCol)
        }
        layoutList[[thisRun]] <- plateDat[order(plateDat$SampleName),c("Row", "Column")]
    }
    
   write.csv(layoutScoreTab, paste(outputFile, plateID, ".PlateDistribution.csv", sep=""), row.names = F, quote=F)
   layoutScoreTab <- apply(layoutScoreTab, 2, function(x){x/mean(x)})
   finalScores <- rowSums(layoutScoreTab)
   whiBest <- which.min(finalScores)
   plateDat <- plateDat[order(plateDat$SampleName),]
   plateDat$Row <- layoutList[[whiBest]]$Row
   plateDat$Column <- layoutList[[whiBest]]$Column
    
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

    if(BatchColumns!="<undefined>"){
        allbatches <- strsplit(BatchColumns, ",")[[1]]
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
            ggsave(paste(outputFile, plateID, ".", colName, ".png", sep=""), height = 5, width=plotWidth)  
        }
    }

    #replace the well row numbers with letters and create the well index column
    plateDat <- plateDat %>% 
        mutate(Row=LETTERS[Row]) %>% 
        mutate(Well=paste(Row, Column, sep=""))

    # bind the plate to the output table
    tabout <- plateDat %>% 
        select(-GroupNumber, -SampleNumber) %>% 
        bind_rows(tabout)
}


# if there is only one plate remove the plate number column. 
tabout <- tabout %>% 
    arrange(PlateNumber, Column, Row)
if(NumberOfPlates==1) dat <- select(dat, -PlateNumber)
#write the final output table
write.csv(tabout, paste(outputFile, ".csv", sep=""), row.names = F, quote=F)


#create a pdf with plots showing the distribution of the different factors across plates
if(NumberOfPlates>1){
    pdf(paste(outputFile, ".Multiplate_Distribution_Check.pdf", sep=""))
    multiplatePlot(tabout, BatchColumns)
    dev.off()
}

############################################################################################
#### NOTES #################################################################################
############################################################################################

# NOTE 1 - Calculating the number of columns and samples on each plate
 # The aim here is to distribute samples as evenly as possible between multiple plates, but also to use the minimum number of columns. 
 # If there are 100 samples, splitting the samples evenly would put 50 samples on each plate. This would use 7 columns on each plate with 12 empty wells (6 on each plate). More efficiently we could put 48 samples on one plate, using 6 full columns, and 52 samples on the other using 7 columns with 4 empty wells.
 # The additional consideration is the Genomics Controls. If the plate is full (12 columns), we need to reserve two wells for Genomics controls.
