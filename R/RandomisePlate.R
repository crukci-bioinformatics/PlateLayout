# This script contains three functions to randomise samples across
# across a 96 well plate given a sample metadata table. The distribution of the samples is optimised using the contents of a single column in the table - the `primaryGroup`.
# Steps:
# 1) randomiseSampleTable - Randomise the order of different factors in the
# "primary group` and replicates within each factor in the metadata table:
#       (1) randomise the order of the sample groups in the table
#       (2) randomise the order of the replicates within each group
# 2) addWells:
#       Assign samples to wells running diagonally across the plate 
#       If the number of columns is even this causes a problem as we cycle back
#       to A1 before completing assignments. We can get around this by adding
#       a dummy column, ordering, and then removing the dummy column
#       We also need to account for empty wells in the final row if the number
#       of samples is not divisible by 8.
# 3) randomiseWells - Randomize the wells by: 
#       (1) randomize the order of the columns except the last one (may be only
#           partially full)
#       (2) randomize the order of the rows
#       (3) move any empty wells in the final column to the bottom

randomiseSampleTable <- function(dat, primaryGroup){
    # (1) randomise the order of the sample groups in the table
    # (2) randomise the order of the replicates within each group
    dat %>% 
        group_by_at(primaryGroup) %>%  
        nest()  %>%  
        ungroup() %>%  
        mutate(Ord.Gp = sample(seq(n()))) %>% 
        unnest(everything()) %>%  
        mutate(Ord.Sam = sample(seq(n()))) %>% 
        arrange(Ord.Gp, Ord.Sam) %>%  
        select(-Ord.Gp, -Ord.Sam)
}


addWells <- function(dat){
    nCols <- ceiling(nrow(dat)/8)
    nCols <- ifelse(nCols%%2==0, nCols+1, nCols) # add dummy column if necessary
    nSams <- nrow(dat)
    tibble(RowNum = rep(1:8, nCols),
           ColNum = rep(seq(nCols), 8)) %>%
           filter(ColNum <= nCols) %>% # remove dummy columns
           rowid_to_column("WellNumber") %>%
           arrange(ColNum, RowNum) %>%  
           slice(n = seq(nSams)) %>%  # remove empty wells in final col
           arrange(WellNumber) %>%  
           select(-WellNumber) %>%
           bind_cols(dat)
}

randomiseWells <- function(dat){ 
    # (1) Randomize columns (except last)
    dat <- dat %>%  
        group_by(ColNum) %>%  
        nest() %>%  
        ungroup() %>%  
        mutate(Column = c(sample(seq(n()-1)), n())) %>%
        unnest(cols = everything()) 
    # (2) Randomize rows
    dat <- dat %>%  
        group_by(RowNum) %>%  
        nest() %>%  
        ungroup() %>%  
        mutate(RowID = sample(seq(n()))) %>%
        arrange(RowID) %>%  
        unnest(cols = everything()) %>%  
        arrange(Column, RowID)
    # (3) Fix empty wells in final column if necessary
    nCols <- ceiling(nrow(dat)/8) # number of columns
    nFinalCol <- nrow(dat)%%8 # number of samples in the final column
    if(nFinalCol>0){
        dat <- dat %>%  
            mutate(RowID = ifelse(Column==nCols, seq(nFinalCol), RowID))
    }
    # (4) Add row letters, add well
    dat %>%  
        mutate(Row = LETTERS[RowID]) %>%  
        select(-ColNum, -RowNum) %>%  
        mutate(Well = str_c(Row, Column))
}


# main function
plateRandomisation <- function(samSht, primaryGroup){
    samSht %>%  
        randomiseSampleTable(primaryGroup = primaryGroup)  %>%  
        addWells() %>%  
        randomiseWells()
}            
