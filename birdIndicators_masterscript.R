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


#---------------#
# Download data #
#---------------#

## Download Trim data, incl. EURING codes, from database
minYear <- 2006
maxYear <- 2023

Trim_data <- downloadData_TRIM(minYear = minYear, maxYear = maxYear,
                               drop_negativeSpp = TRUE)

# NOTE: What gets downloaded from the database here are Trim RESULTS. 
# How were/are these generated? Does this happen internally in the database
# or is it done manually?


#-----------------------#
# PECBMS Analysis setup #
#-----------------------#

## Get species lists
sppLists <- makeSpeciesLists(Trim_data = Trim_data)
Spp_selection <- sppLists$sppData

## Set directory
folder <- "PECBMS_Files"

## Write PECBMS arguments input files for each species
argument_file <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection,
                                                  folderPath = folder)

## Subset data to contain only relevant species
PECBMS_data <- makeInputData_PECBMS(Trim_data = Trim_data,
                                    Spp_selection = Spp_selection,
                                    convertNA = TRUE, 
                                    save_allSppData = TRUE, returnData = TRUE)

#------------------#
# PECBMS Trim runs #
#------------------#

## Run analyses using the PECBMS Rtrim shell
runRtrimShell_PECBMS(folder = folder)


#-----------------------------#
# Process PECBMS Trim results #
#-----------------------------#

## Process results
trimResults_PECBMS <- processRtrimOutput_PECBMS(folder = folder)


#--------------------------------#
# Post-processing: collect files #
#--------------------------------#

## Set name of folder in which to collect files
subFolderName <- "Species_files"

## Collect (and rename) species and summary files
collectSpeciesFiles_PECBMS(folder = folder, 
                           subFolderName = subFolderName)

## Retrieve equivalent file from predecessor monitoring programme
legacyFile_folder <- "P:/41201612_naturindeks_2021_2023_database_og_innsynslosning/Hekkefugl_Dataflyt/LegacyFiles_PECBMS_Trim"
collectSpeciesFiles_Legacy(origin_folder = legacyFile_folder,
                           target_folder = paste0(folder, "/", subFolderName))

#------------------------#
# PECBMS RSWAN execution #
#------------------------#

## Set working folder
working_folder = paste0(folder, "/", subFolderName)

## Write/load schedule table
writeSchedule_SWAN(working_folder = working_folder, 
                   general_folder = "data",
                   MSI_speciesList = sppLists$sppLists$MSI,
                   loadSchedule = FALSE)

