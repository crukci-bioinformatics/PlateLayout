% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainFunctions.R
\name{randomizeSinglePlate}
\alias{randomizeSinglePlate}
\title{Randomize Single Plate}
\usage{
randomizeSinglePlate(designTable, outputFile = NULL,
  primaryGroup = "SampleGroup", batchColumns = NULL, nIter = 10000,
  nCores = 4)
}
\arguments{
\item{designTable}{The meta data table - a data.frame/tibble.}

\item{outputFile}{The root name for the output files if required.}

\item{primaryGroup}{The name of a column in the metadata table to use for
optimising the distribution of samples across the plate}

\item{batchColumns}{A character vector wit the names of additional columns
in the metadata table to use in assessing the distribution of the samples
across the plate.}

\item{nIter}{The number of random layouts to generate and assess.}

\item{nCores}{The number of cores to use. Set to 1 to run serial.}
}
\value{
If no \code{outputFile} is provided, a tibble containing the original
metadata and the well assignments for each sample. If an \code{outputFile} is
provided, nothing - the function writes out the plate layout to a tsv file 
and a plot for each of `primaryGroup` and `batchColumns` to png files.
}
\description{
Randomize up to 96 samples across a single 96 well plate.
}
\details{
The metadata (\code{designTable}) should be a data.frame/tibble. The 
column names can be anything. By default the function expects one column 
named \code{SampleGroup} that it will use to distribute the samples across
the plate. This can be set to any other column using the \code{primaryGroup}
argument.

The columns specified in \code{batchColumns} are used in addition
to \code{primaryGroup} to assess which of the \code{nIter} layouts is 
optimal in terms of the distribution of samples. This is assessed based on
the sum of Euclidean distances between replicate samples within each column.
}
\examples{
library(readr)
designSheet <- system.file("extdata", "metadata_12x3.tsv",
                           package = "PlateLayout")
bColumns <- c("ExtractionInformation", "PassageNumber") 
designTable <- read_tsv(designSheet)
# Create an R object in the session
plateLayout <- randomizeSinglePlate(designTable,  batchColumns = bColumns,
                     nIter = 100) 

# Or directly output a table and plots
outputFile <- "Test"
randomizeSinglePlate(designTable, 
                     outputFile = outputFile, 
                     batchColumns = bColumns,
                     nIter = 100) 
}
