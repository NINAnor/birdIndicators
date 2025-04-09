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

if(!dir.exists("NI_regions_files")){dir.create("NI_regions_files")}

folderN <- "NI_regions_files/Species_filesN"
folderS <- "NI_regions_files/Species_filesS"
folderS_only <- "NI_regions_files/Species_filesS_only"
folderGrouse <- "NI_regions_files/Species_files_Grouse"
folderHele <- "NI_regions_files/Species_files_HeleNorge"

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


#--------------------------------------------------#
# Add the regional differences - subset the routes #
#--------------------------------------------------#

## Add a column with the specifications for N and S data
Trim_data <- Trim_data %>%
  mutate(Regions = case_when(Site %in% c(101:1736) ~ 'South',
                             TRUE ~'North'))



#****************************************************************************
#
# PECBMS MAIN TRIM & SWAN RUNS
#
#****************************************************************************

#---------------------------------#
# Main run: PECBMS Analysis setup #
#---------------------------------#

## Get species lists (and save as file)
sppLists <- makeSpeciesLists(Trim_data = Trim_data)
saveRDS(sppLists, file = "data/sppLists.rds")

## Make species selection for main PECBMS runs
Spp_selection <- sppLists$sppData %>%
  dplyr::filter(WEBSITE)

## Write PECBMS arguments input files for each species
argument_file <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection,
                                                  folderPath = folder)

## Subset data to contain only relevant species
PECBMS_data <- makeInputData_PECBMS(Trim_data = Trim_data,
                                    Spp_selection = Spp_selection,
                                    convertNA = TRUE, 
                                    save_allSppData = TRUE, returnData = TRUE)

#----------------------------#
# Main run: PECBMS Trim runs #
#----------------------------#

## Run analyses using the PECBMS Rtrim shell
runRtrimShell_PECBMS(folder = folder)


#---------------------------------------#
# Main run: Process PECBMS Trim results #
#---------------------------------------#

## Process results
trimResults_PECBMS <- processRtrimOutput_PECBMS(folder = folder)


#------------------------------------------#
# Main run: Post-processing: collect files #
#------------------------------------------#

## Collect (and rename) species and summary files
collectSpeciesFiles_PECBMS(folder = folder, 
                           subFolderName = subFolderName)

## Retrieve equivalent file from predecessor monitoring programme (legacy files)
collectSpeciesFiles_Legacy(origin_folder = legacyFile_folder,
                           target_folder = paste0(folder, "/", subFolderName))

#--------------------------#
# PECBMS RSWAN preparation #
#--------------------------#

## Make sub-directory if not already present
if(!file.exists(output_folder2)){ 
  dir.create(output_folder2)
}

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



#****************************************************************************
#
# MULTISPECIES INDICES
#
#****************************************************************************

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


#****************************************************************************
#
# NATURE INDEX
#
#****************************************************************************

#-------------------------------------#
# Nature Index: PECBMS Analysis setup #
#-------------------------------------#

## For NI: Select species where indexes have to be calculated with data split into N and S datasets
Spp_selection_N_S <- sppLists$sppData %>%
  filter(dataUse_direct_NN & dataUse_direct_SN)

## For NI: Species where only data from S are used
Spp_selection_S_only <- sppLists$sppData %>%
  filter(dataUse_direct_SN & !dataUse_direct_NN)

## For NI: Species where data from the whole country are used.
Spp_selection_HeleNorge <- sppLists$sppData %>%
  filter(dataUse_direct_N | dataUse_expert_N)

## For NI: Gallinaceous - Data split into four regions
Spp_selection_4 <- sppLists$sppData %>%
  filter(dataType == "TBD")

## For NI: Subset TRIM data to calculate indexes for N and S
Trim_data_N <- Trim_data %>%
  filter(Spp_name %in% Spp_selection_N_S$Species & Regions == 'North') 

Trim_data_S <- Trim_data %>%
  filter(Spp_name %in% Spp_selection_N_S$Species & Regions == 'South')

Trim_data_S_only <- Trim_data %>%
  filter(Spp_name %in% Spp_selection_S_only$Species & Regions == 'South')

Trim_data_HeleNorge <- Trim_data %>%
  filter(Spp_name %in% Spp_selection_HeleNorge$Species)

Trim_data_grouses <- Trim_data %>%
  filter(Spp_name %in% Spp_selection_4$Species)


## Write PECBMS arguments input files for each species by regions of interest to NI
argument_file_N <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection_N_S,
                                                    folderPath = folderN)

argument_file_S <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection_N_S,
                                                    folderPath = folderS)

argument_file_S_only <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection_S_only,
                                                         folderPath = folderS_only)

argument_file_HeleNorge <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection_HeleNorge,
                                                             folderPath = folderHele)

argument_file_Grouse <- setupInputFiles_PECBMS_trimShell(Spp_selection = Spp_selection_4,
                                                         folderPath = folderGrouse)

## Subset data to contain only relevant species and regions
NI_data_N <- makeInputData_NI_N(Trim_data = Trim_data_N,
                                Spp_selection = Spp_selection_N_S,
                                convertNA = TRUE, 
                                save_allSppData = TRUE, returnData = TRUE)

NI_data_S <- makeInputData_NI_S(Trim_data = Trim_data_S,
                                Spp_selection = Spp_selection_N_S,
                                convertNA = TRUE, 
                                save_allSppData = TRUE, returnData = TRUE)

NI_data_S_only <- makeInputData_NI_S_only(Trim_data = Trim_data_S_only,
                                          Spp_selection = Spp_selection_S_only,
                                          convertNA = TRUE, 
                                          save_allSppData = TRUE, returnData = TRUE)

NI_data_HeleNorge <- makeInputData_NI_HeleNorge(Trim_data = Trim_data_HeleNorge,
                                                Spp_selection = Spp_selection_HeleNorge,
                                                convertNA = TRUE, 
                                                save_allSppData = TRUE, returnData = TRUE)

NI_data_Grouse <- makeInputData_NI_Grouse(Trim_data = Trim_data_grouses,
                                          Spp_selection = Spp_selection_4,
                                          convertNA = TRUE, 
                                          save_allSppData = TRUE, returnData = TRUE)


#-------------------------------------------#
# Nature Index: PECBMS Trim runs per region #
#-------------------------------------------#

## Run analyses using the PECBMS Rtrim shell
runRtrimShell_PECBMS(folder = folderN)

runRtrimShell_PECBMS(folder = folderS)

runRtrimShell_PECBMS(folder = folderS_only) 

runRtrimShell_PECBMS(folder = folderHele)

runRtrimShell_PECBMS(folder = folderGrouse) 


#-------------------------------------------------#
# Nature Index: Process NI Trim results by region #
#-------------------------------------------------#

## Process results
trimResults_NI_N <- processRtrimOutput_PECBMS(folder = folderN)

trimResults_NI_S <- processRtrimOutput_PECBMS(folder = folderS)

trimResults_NI_S_only <- processRtrimOutput_PECBMS(folder = folderS_only)

trimResults_NI_HeleNorge <- processRtrimOutput_PECBMS(folder = folderHele)

trimResults_NI_Grouse <- processRtrimOutput_PECBMS(folder = folderGrouse)


#----------------------------------------------#
# Nature Index: Post-processing: collect files #
#----------------------------------------------#

## Collect (and rename) species and summary files
collectSpeciesFiles_PECBMS(folder = folderN, 
                           subFolderName = subFolderName)

collectSpeciesFiles_PECBMS(folder = folderS, 
                           subFolderName = subFolderName)

collectSpeciesFiles_PECBMS(folder = folderS_only, 
                           subFolderName = subFolderName)

collectSpeciesFiles_PECBMS(folder = folderHele, 
                           subFolderName = subFolderName)

collectSpeciesFiles_PECBMS(folder = folderGrouse, 
                           subFolderName = subFolderName)

#------------------------------------#
# Nature index indicator calculation #
#------------------------------------#

## Set parameters
NI_years <- c(2000, 2010, 2014, 2019, 2024)
refAnchorYear <- 2010 # Reference anchor year
refAnchorYear_alt <- 2019 # Alternative reference anchor year (for indicators without data in 2010)
nsim <- 10000 # Number of simulations for averaging


## Register NI database credentials & request token
UserName_NIdb <- rstudioapi::askForPassword("NI database username") # = NINA email address
Password_NIdb <- rstudioapi::askForPassword("NI database password")

NIcalc::getToken(username = UserName_NIdb,  
                 password = Password_NIdb,
                 url = "https://www8.nina.no/NaturindeksNiCalc")

## Re-load species list if not present
if(!exists("sppLists")){
  sppLists <- readRDS("data/sppLists.rds")
}
Spp_selection <- sppLists$sppData

## Categorise species flagged for NI according to their updating approach
expertJudge <- Spp_selection$Species[which(Spp_selection$dataUse_expert_N)]

otherUse <- Spp_selection$Species[which(Spp_selection$dataType == "TBD")]

## List species relevant for upload to NI database (directly or via expert assessment)
sppNI <- listSpecies_NI(Spp_selection = Spp_selection,
                        Spp_exclude = otherUse)


## Download old indicator data from NI database and store as a list
OldIndicator_data <- list()
for(i in 1:nrow(sppNI)){
  OldIndicator_data[[i]] <- NIcalc::getIndicatorValues(indicatorID = sppNI$indicatorId[i])
}
names(OldIndicator_data) <- sppNI$indicatorId

## List potentially relevant input file folders
inputFile_folders <- data.frame(
  ID = c("PECBMS",
         "NI_Norge",
         "NI_NNorge",
         "NI_SNorge",
         "NI_onlySNorge"),
  path = c(working_folder_rel,
           paste0(folderHele, "/Species_files"),
           paste0(folderN, "/Species_files"),
           paste0(folderS, "/Species_files"),
           paste0(folderS_only, "/Species_files"))
)

## Determine "updating plan" (which Trim file for which NI area) for all species
updatePlans <- makeUpdatePlan_NI(sppNI = sppNI,
                                 OldIndicator_data = OldIndicator_data)


## Assemble and average TRIM data & use it to calculate NI indicator data for each species-area
UpdatedIndicator_data <- prepareIndicatorData_NI(sppNI = sppNI,
                                                 OldIndicator_data = OldIndicator_data,
                                                 updatePlans = updatePlans,
                                                 inputFile_folders = inputFile_folders,
                                                 use_combTS = TRUE,
                                                 NI_years = NI_years, 
                                                 refAnchorYear = refAnchorYear, 
                                                 refAnchorYear_alt = refAnchorYear_alt,
                                                 nsim = nsim)

## Write updated indicator data to csv for review
writeIndicatorData_forReview(UpdatedIndicator_data = UpdatedIndicator_data,
                             dir = "NI_indicatorData_forReview",
                             sppNI = sppNI)

## Upload updated indicator data for select species to NI database

indicatorId_upload <- c(278, 289, 301, 305, 308, 313, 314, 317, 318, 325, 326, 
                        328, 330, 335, 389, 392:403, 405:410, 412, 414, 415 )


for(i in 1:length(UpdatedIndicator_data)){
  
  if(names(UpdatedIndicator_data)[i] %in% indicatorId_upload){
    message(paste0("Upload data for indicator: ", names(UpdatedIndicator_data)[i]))
    NIcalc::writeIndicatorValues(UpdatedIndicator_data[[i]])
  }else{
    message(paste0("Skipping indicator: ", names(UpdatedIndicator_data)[i]))
  }

}  
