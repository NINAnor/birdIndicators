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
library(readxl)
library(data.table)
library(NIcalc) # For Nature Index only; manual download from GitHub

#----------------#
# Workflow setup #
#----------------#

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

## Set relative directory/file paths

# Storage of raw PECBMS Trim outputs
folder <- "PECBMS_Files" 

# Storage of processed PECBMS Trim outputs
subFolderName <- "Species_files" 

# Relative general and working folders for RSWAN
general_folder_rel <- "data"
working_folder_rel <- paste0(folder, "/", subFolderName)

# Storage of MSI results
MSI_results_folder <- "MSI_Results"


## Set absolute directory/file paths

# NINAs local storage of legacy (1995 - 2008) Trim output files
legacyFile_folder <- "P:/41201612_naturindeks_2021_2023_database_og_innsynslosning/Hekkefugl_Dataflyt/LegacyFiles_PECBMS_Trim"

# Absolute data and output paths for RSWAN 
# (PECBMS' RSWAN scripts are not compatible with relative paths as they 
# repeatedly use setwd() to toggle between folders)
general_folder <- paste0(getwd(), "/", general_folder_rel)
working_folder <- paste0(getwd(), "/", working_folder_rel)

output_folder <- paste0(working_folder, "/output/") 
output_folder2 <- paste0(working_folder, "/output") 

if(!file.exists(output_folder2)){ 
  dir.create(output_folder2)
}
# NOTE: The RSWAN scripts require that absolute folder paths are assigned to 
# objects strictly named "general_folder", "working_folder", and "output_folder"
# because some of the RSWAN helper functions use setwd() calls that refer to 
# globally defined absolute paths. 



#---------------#
# Download data #
#---------------#

## Download Trim data, incl. EURING codes, from database
minYear <- 2006
maxYear <- 2024

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

## Collect (and rename) species and summary files
collectSpeciesFiles_PECBMS(folder = folder, 
                           subFolderName = subFolderName)

## Retrieve equivalent file from predecessor monitoring programme (legacy files)
collectSpeciesFiles_Legacy(origin_folder = legacyFile_folder,
                           target_folder = paste0(folder, "/", subFolderName))

#--------------------------#
# PECBMS RSWAN preparation #
#--------------------------#

## Write/load schedule table
writeSchedule_SWAN(working_folder = working_folder_rel, 
                   general_folder = general_folder_rel,
                   MSI_speciesList = sppLists$sppLists$MSI,
                   loadSchedule = FALSE)

## Truncate first survey year where necessary
correctFirstSurveyYear_SWAN(general_folder = general_folder_rel, 
                            working_folder = working_folder_rel,
                            maxYear = maxYear)

#------------------------#
# PECBMS RSWAN execution #
#------------------------#

## Run RSWAN
combineTimeSeries_SWAN(general_folder_abs = general_folder,
                       working_folder_abs = working_folder)



#--------------------------------------#
# Multispecies Index (MSI) calculation #
#--------------------------------------#

## Determine data range, base (= reference) year, and changepoint year
useCombTS <- TRUE
baseYear <- 2000 # Previously 1996
changepointYear <- 2008

## Attempt to calculate long-term MSIs for four ecosystems
IndexNames <- c("FarmlandBirds", "ForestBirds", "MountainBirds", "WetlandBirds")
Index_sppLists <- list(sppLists$sppLists$MSI_farmland,
                       sppLists$sppLists$MSI_forest,
                       sppLists$sppLists$MSI_mountain,
                       sppLists$sppLists$MSI_wetlands)

MSI_longTerm <- data.frame()
for(i in 1:2){
  message(paste0(crayon::bold("Calculating ", IndexNames[i], " Index...")))
  results <- calculateIndex_MultiSpecies(working_folder = working_folder_rel, 
                                         Spp_subset = Index_sppLists[[i]], 
                                         IndexName = IndexNames[i], 
                                         results_folder = MSI_results_folder,
                                         useCombTS = useCombTS,
                                         baseYear = baseYear,
                                         changepointYear = changepointYear)
  MSI_longTerm <- rbind(MSI_longTerm, results)
  message("")
}

## Calculate long-term and mid-term MSIs for all ecosystems
useCombTS <- c(TRUE, TRUE, FALSE, FALSE)
baseYear <- 2008

MSI_all <- data.frame()
for(i in 1:length(IndexNames)){
  message(paste0(crayon::bold("Calculating ", IndexNames[i], " Index...")))
  results <- calculateIndex_MultiSpecies(working_folder = working_folder_rel, 
                                            Spp_subset = Index_sppLists[[i]], 
                                            IndexName = IndexNames[i], 
                                            results_folder = MSI_results_folder,
                                            useCombTS = useCombTS[i],
                                            baseYear = baseYear,
                                            changepointYear = changepointYear)
  MSI_all <- rbind(MSI_all, results)
  message("")
}


## Save results as .rds
MSI_results <- list(MSI_baseline2000 = MSI_longTerm,
                    MSI_baseline2008 = MSI_all)

saveRDS(MSI_results, file = paste0(MSI_results_folder, "/MSI_results.rds"))


#----------------------------------#
# Plot Multispecies Indices (MSIs) #
#----------------------------------#

## Plot farmland and forest MSIs with baseline 1996
plotTimeSeries_MSI(MSI = MSI_results$MSI_baseline2000,
                   results_folder = MSI_results_folder,
                   plot_name = "MSI_base2000_AllEcosystems", 
                   displayPlots = TRUE,
                   savePDF = TRUE)


## Plot all MSIs with baseline 2008
plotTimeSeries_MSI(MSI = MSI_results$MSI_baseline2008,
        results_folder = MSI_results_folder,
        plot_name = "MSI_base2008_AllEcosystems",
        displayPlots = TRUE,
        savePDF = TRUE)


#------------------------------------#
# Nature index indicator calculation #
#------------------------------------#

## Set parameters
NI_years <- c(2000, 2010:2014, 2019, 2024)
refAnchorYear <- 2010 # Reference anchor year
nsim <- 10000 # Number of simulations for averaging


## Register NI database credentials & request token
UserName_NIdb <- rstudioapi::askForPassword("NI database username") # = NINA email address
Password_NIdb <- rstudioapi::askForPassword("NI database password")

NIcalc::getToken(username = UserName_NIdb,  
                 password = Password_NIdb,
                 url = "https://www8.nina.no/NaturindeksNiCalc")

## Categorise species flagged for NI according to their updating approach
expertJudge <- c("Acrocephalus schoenobaenus",
                 "Anthus petrosus",
                 "Aythya marila",
                 "Calidris alpina",
                 "Charadrius morinellus",
                 "Clangula hyemalis",
                 "Eremophila alpestris",
                 "Gallinago media",
                 "Melanitta fusca",
                 "Melanitta nigra",
                 "Phalaropus lobatus",
                 "Plectrophenax nivalis" 
                 )

otherUse <- c("Lagopus lagopus", 
              "Lagopus muta", 
              "Lyrurus tetrix", 
              "Tetrao urogallus")

directNI <- listSpecies_NI(Spp_selection = Spp_selection,
                           Spp_exclude = c(expertJudge, otherUse))



## Assemble and average TRIM data & use it to calculate NI indicator data for each species-area
IndData <- prepareIndicatorData_NI(directNI = directNI,
                                   working_folder = working_folder)


#TODO: Once we have scripted the subpopulations for TRIM analyses, we have to
# update the above function to correctly deal with the case of several areas 
# with distinct values)

#TODO: Once we add the area-specific data, we have to go over 
# IndData$UpdatedIndicator_data and "merge" the data for different areas of the 
# same species.

## Write updated indicator data to csv for review
writeIndicatorData_forReview(UpdatedIndicator_data = IndData$UpdatedIndicator_data,
                             dir = "NI_indicatorData_forReview")

