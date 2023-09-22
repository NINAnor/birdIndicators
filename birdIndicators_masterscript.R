library(tidyverse)
library(rtrim)
library(data.table)
library(RPostgres)
library(odbc)
library(RPostgreSQL)
library(DBI)
library(rpostgis)
library(sf)
library(lubridate)
library(xtable)
library(rtrim)

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

#-----------------------#
# PECBMS Analysis setup #
#-----------------------#

## Get the selected 71 species and years to be included in the report
file_sppInfo_PECBMS <- "data/ReportSpeciesPECBMS2020.rds"
file_sppList_PECBMS <- "data/PECBMS_species_list_2022.rds"

Spp_selection <- listSpecies_PECBMS(file_sppInfo_PECBMS = file_sppInfo_PECBMS, 
                                    file_sppList_PECBMS = file_sppList_PECBMS)


## Write PECBMS arguments input files for each species
argument_file <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection,
                                                  folderPath = "PECBMS_Files")


#-----------------#
# TOV-E data prep #
#-----------------#

## Download Trim data, incl. EURING codes, from database
minYear <- 2006
maxYear <- 2023

Trim_data <- downloadData_TRIM(minYear = minYear, maxYear = maxYear,
                               drop_negativeSpp = TRUE)

# NOTE: What gets downloaded from the database here are Trim RESULTS. 
# How were/are these generated? Does this happen internally in the database
# or is it done manually?

## Subset data to contain only relevant species
PECBMS_data <- makeInputData_PECBMS(Trim_data = Trim_data,
                                    Spp_selection = Spp_selection,
                                    convertNA = TRUE, 
                                    save_allSppData = TRUE, returnData = TRUE)



#------------------#
# PECBMS Trim runs #
#------------------#

## Set directory
folder <- "PECBMS_Files"

## Run analyses using the PECBMS Rtrim shell
runRtrimShell_PECBMS(folder = folder)


#-----------------------------#
# Process PECBMS Trim results #
#-----------------------------#

## Process results
trimResults_PECBMS <- processRtrimOutput_PECBMS(folder = folder)




subFolderName <- "Species_files"