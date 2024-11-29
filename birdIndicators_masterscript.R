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
baseYear <- 1996
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
MSI_results <- list(MSI_baseline1996 = MSI_longTerm,
                    MSI_baseline2008 = MSI_all)

saveRDS(MSI_results, file = paste0(MSI_results_folder, "/MSI_results.rds"))


#----------------------------------#
# Plot Multispecies Indices (MSIs) #
#----------------------------------#

## Plot farmland and forest MSIs with baseline 1996
plotTimeSeries_MSI(MSI = MSI_results$MSI_baseline1996,
                   results_folder = MSI_results_folder,
                   plot_name = "MSI_base1996_AllEcosystems", 
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

directNI <- Spp_selection %>%
  dplyr::filter(NI & !(Species %in% c(expertJudge, otherUse)))


## Add Norwegian names for Nature index species
NorwegianNames <- read.csv("data/NorwegianSppNames.csv")
directNI <- directNI %>%
  dplyr::left_join(NorwegianNames, by = c("EURINGCode", "Species"))


# Check for missing species
if(any(is.na(directNI$Species_Norwegian_NI))){
  NAs <- subset(directNI, is.na(Species_Norwegian_NI))
  message(paste0("Norwegian species names are missing for the following ", nrow(NAs), " species:"))
  print(NAs$Species)
  
  message("If these are new species to be included in NI, please add their Norwegian name to the file data/NorwegianSppNames.csv")
}

if(!all(NorwegianNames$Species %in% directNI$Species)){
  NAs <- subset(NorwegianNames, !(Species %in% directNI$Species))
  message(paste0("The list of Norwegian names includes", nrow(NAs), " species that are part of Nature Index but for which there seems to be no entry in directNI:"))
  print(NAs$Species)
  
  message("Please ensure that all necessary species are included in directNI (check dependencies of function makeSpeciesLists).")
  
}


## Register NI database credentials & request token
UserName_NIdb <- rstudioapi::askForPassword("NI database username") # = NINA email address
Password_NIdb <- rstudioapi::askForPassword("NI database password")

NIcalc::getToken(username = UserName_NIdb,  
                 password = Password_NIdb,
                 url = "https://www8.nina.no/NaturindeksNiCalc")


## Retrieve NI indicator IDs for all relevant species
myIndicators <- NIcalc::getIndicators() %>%
  dplyr::rename(Species_Norwegian_NI = name)

directNI <- directNI %>%
  dplyr::left_join(myIndicators, by = "Species_Norwegian_NI")

# Check that all relevant species have an id in the Nature Index database
if(any(is.na(directNI$id))){
  NAs <- subset(directNI, is.na(id))
  message(paste0("Indicator ID is missing for the following ", nrow(NAs), " species:"))
  print(NAs$Species)
  
  message("If these are new species to be included in NI, you will have to request their addition to the Nature Index database.")
  message("If these species already exist in the NI database, you may be missing indicator access.")
}


## Set up lists for storing data
TrimIndex_data <- list()
OldIndicator_data <- list()
UpdatedIndicator_data <- list()


## For each species:
for(i in 1:nrow(directNI)){
  
  message(paste0("Preparing data for indicator: ", directNI$Species_Norwegian_NI[i]))
  
  # Download old indicator data from NI database
  message("Retrieving data from NI database...")
  OldIndicator_data[[i]] <- NIcalc::getIndicatorValues(indicatorID = directNI$id[i])

  # Check areas used in NI database
  areas <- unique(OldIndicator_data[[i]]$indicatorValues$areaName)
  message("NI database contains the following ", length(areas), " area(s):")
  print(areas)
  
  #TODO: Once we have scripted the subpopulations for TRIM analyses, we have to
  # update this part to correctly deal with the different cases (one area, 
  # several areas but same values, several areas with distinct values)
  
  if(length(areas) > 1){
    message("Indicator is defined for mulitple areas. Functionality for this is not implemented yet, but will be added soon. For now, this indicator will be skipped.")
    next()
  }
  
  # Check for availability of (combined) TRIM index data
  spp_files <- list.files(working_folder)[which(grepl(directNI$EURINGCode[i], list.files(working_folder)))]
  combTS <- ifelse(any(grepl("COMB", spp_files)), TRUE, FALSE)
  
  # Load relevant TRIM index data 
  message("Loading TRIM index data...")
  if(combTS){
    index_TS_file <- spp_files[which(grepl("indices_TT", spp_files) & grepl("COMB", spp_files))]
  }else{
    index_TS_file <- spp_files[which(grepl("indices_TT", spp_files))]
  }
  
  TrimIndex_data[[i]] <- read.csv(paste0(working_folder, "/", index_TS_file), sep = ";")
  
  # Calculate 3-year averages for TRIM index data
  message("Calculating 3-year averages for TRIM data...")
  TrimIndex_averages <- data.frame()
  
  for(t in 3:nrow(TrimIndex_data[[i]])){
    idx <- matrix(NA, nrow = nsim, ncol = 3)
    
    for(s in 1:3){
      idx[, s] <- truncnorm::rtruncnorm(n = nsim, 
                                        mean = TrimIndex_data[[i]]$Index_model[t-(s-1)],
                                        sd = TrimIndex_data[[i]]$Index_model_SE[t-(s-1)],
                                        a = 0)
    }
    
    avg_3yr <- rowMeans(idx)
    
    avg_data <- data.frame(Year = TrimIndex_data[[i]]$Year[t],
                           mean = mean(avg_3yr),
                           sd = sd(avg_3yr),
                           lowerQuart = unname(quantile(avg_3yr, probs = 0.25)),
                           upperQuart = unname(quantile(avg_3yr, probs = 0.75)))
    
    TrimIndex_averages <- rbind(TrimIndex_averages, avg_data)
  }
  
  
  # Extract reference proportion for the reference anchor year from old NI data
  refProp_orig <- subset(OldIndicator_data[[i]]$indicatorValues, yearName == refAnchorYear)$verdi/100
  
  # Calculate and add reference value
  message("Recalibrate reference value...")
  ref <- subset(TrimIndex_averages, Year == refAnchorYear) %>%
    dplyr::mutate(mean = mean/refProp_orig,
                  sd = sd/refProp_orig,
                  lowerQuart = lowerQuart/refProp_orig,
                  upperQuart = upperQuart/refProp_orig)
  ref[,"Year"] <- "Referanseverdi"
  
  TrimIndex_averages <- rbind(TrimIndex_averages, ref)
  
  # Do pre-scaling and drop values not relevant for upload to NI database
  message("Scale and format updated indicator data...")
  Indicator_prescaled <- TrimIndex_averages %>%
    dplyr::mutate(mean = 100*mean/ref$mean,
                  sd = 100*sd/ref$mean,
                  lowerQuart = 100*lowerQuart/ref$mean,
                  upperQuart = 100*upperQuart/ref$mean) %>%
    dplyr::filter(Year %in% c("Referanseverdi", NI_years)) %>%
    dplyr::rename(yearName = Year)
  
  # Write updated indicator data in upload format
  Indicator_upload <- OldIndicator_data[[i]]$indicatorValues %>%
    dplyr::left_join(Indicator_prescaled, by = "yearName") %>%
    dplyr::mutate(update = ifelse(is.na(mean), FALSE, TRUE)) %>%
    dplyr::mutate(
      
      # (Over)write indicator values
      verdi = ifelse(update, mean, verdi),
      nedre_Kvartil = ifelse(update, lowerQuart, nedre_Kvartil),
      ovre_Kvartil = ifelse(update, upperQuart, ovre_Kvartil),
      
      # Set data type for updated values
      datatypeId = ifelse(update & yearId != "Referanseverdi", 3, datatypeId),
      datatypeName = ifelse(update & yearId != "Referanseverdi", "Beregnet fra modeller", datatypeId),
      
      # Remove no longer valid information on distributions
      distributionName = ifelse(update, NA, distributionName),
      distributionId = ifelse(update, NA, distributionId),
      distParam1 = ifelse(update, NA, distParam1),
      distParam2 = ifelse(update, NA, distParam2)) %>%
    
    # Delete auxiliary columns
    dplyr::select(-mean, -sd, -lowerQuart, -upperQuart, -update)
  
  UpdatedIndicator_data[[i]] <- OldIndicator_data[[i]]
  UpdatedIndicator_data[[i]]$indicatorValues <- Indicator_upload
  
  # Check for any legacy custom distribution information
  if(length(UpdatedIndicator_data[[i]]$customDistributions) > 0){
    warning(paste0("Legacy custom distribution information is present for species ", directNI$Species_Norwegian_NI[i], ". Please review and delete information."))
  }
}

## Write updated indicator data to csv for review
if(!dir.exists(file.path("NI_indicatorData_forReview"))){
  dir.create(file.path("NI_indicatorData_forReview"))
}

for(i in 1:length(UpdatedIndicator_data)){
  
  if(!is.null(UpdatedIndicator_data[[i]])){
    indID <- UpdatedIndicator_data[[i]]$indicatorValue$indicatorId[1]
  
    readr::write_excel_csv(UpdatedIndicator_data[[i]]$indicatorValues,
                           file = paste0("NI_IndicatorData_forReview/IndID_", indID, ".csv"))
  }
}
