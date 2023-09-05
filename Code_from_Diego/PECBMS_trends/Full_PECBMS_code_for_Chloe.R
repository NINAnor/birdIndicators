################################################################################
############ FULL SCRIPT TO CREATE PECBMS FILES TO CHLO?? #######################
################################################################################
### Version May 2022 ###


require(tidyverse)
require(rtrim)
require(data.table)


folder <- "C:/Users/diego.pavon-jordan/OneDrive - NINA/Documents/PECBMS/2023_Analysis"  
setwd(folder)


## Get the selected 71 species and years to be included in the report
Species_info <- read.csv2('ReportSpeciesPECBMS2020.csv')

Spp_selection <- read.csv2('PECBMS_species_list_2022.csv') %>%
  left_join(Species_info) %>%
  select(EURINGCode, Year_First, TCF_Basetime, TCF_Slope_From) %>%
  as_tibble()%>%
  mutate(Stratum = 63)
# 71 species - Ok!

## As per JAK email -> We use Year_First from ReportSpeciesPECBMS2020 and
## set TCF_Slope_From as the same year as Year_First

Spp_selection3 <- Spp_selection %>%
  mutate(TCF_Slope_From_NEW = Year_First)
# all species have basetime = 2008. Now, Year_First and TCF_Slope_From_NEW
# have the same value


## so I split the data set in those species with trend
## starting in 2007 (and 2006) and those starting in 2008 (based on TCF_Slope_From_NEW)


## Group species by starting year of time series

## 2007 (and 2006)
Spp_selection_year07 <- Spp_selection3 %>%
  filter(TCF_Slope_From_NEW == 2006 | TCF_Slope_From_NEW == 2007)
# 54 species

euring_codes07 <- Spp_selection_year07$EURINGCode

unique(Spp_selection_year07$Year_First)
# 2006 and 2007


## 2008
Spp_selection_year08 <- Spp_selection3 %>%
  filter(TCF_Slope_From_NEW == 2008)
# 17 species

euring_codes08 <- Spp_selection_year08$EURINGCode

unique(Spp_selection_year08$Year_First)
# 2008


## Get the file names as required by PECBMS
File_names07 <- NULL
File_names08 <- NULL

for(i in 1:length(euring_codes07)){
  File_names07[i] <- paste0('BMP_', euring_codes07[i],'_1_63')
}

for(i in 1:length(euring_codes08)){
  File_names08[i] <- paste0('BMP_', euring_codes08[i],'_1_63')
}

File_names07
File_names08



## Saving the species-specific input files as per rtrim-shell format request
## For trends starting in 2007

File <- NA
Base_year_first_year <- 2008
Base_year_last_year <- 2008
Changepoints <- 'all'
Serial_correlation <- 'TRUE'
Overdispersion <- 'TRUE'
Presence_weights <- 'FALSE'
Presence_monthfactors <- 'FALSE'
Year_from <- 2007
Save_fitted_values <- 'TRUE'

argument_file <- data.frame(File, Base_year_first_year, Base_year_last_year,
                            Changepoints, Serial_correlation, Overdispersion,
                            Presence_weights, Presence_monthfactors, Year_from,
                            Save_fitted_values)

for(i in 1:length(File_names07)){
  
  argument_file[i,1] <- File_names07[i]
  argument_file[i,2] <- 2008
  argument_file[i,3] <- 2008
  argument_file[i,4] <- 'all'
  argument_file[i,5] <- 'TRUE'
  argument_file[i,6] <- 'TRUE'
  argument_file[i,7] <- 'FALSE'
  argument_file[i,8] <- 'FALSE'
  argument_file[i,9] <- 2007
  argument_file[i,10] <- 'TRUE'
  
  write.csv2(argument_file[i, ], paste0('BMP_', euring_codes07[i],'_1_63_arg_input_stratum.csv'), row.names = FALSE)
}



## And for 2008

File <- NA
Base_year_first_year2 <- 2008
Base_year_last_year2 <- 2008
Changepoints <- 'all'
Serial_correlation <- 'TRUE'
Overdispersion <- 'TRUE'
Presence_weights <- 'FALSE'
Presence_monthfactors <- 'FALSE'
Year_from <- 2008
Save_fitted_values <- 'TRUE'

argument_file2 <- data.frame(File, Base_year_first_year2, Base_year_last_year2,
                             Changepoints, Serial_correlation, Overdispersion,
                             Presence_weights, Presence_monthfactors, Year_from,
                             Save_fitted_values)


for(i in 1:length(File_names08)){
  
  argument_file2[i,1] <- File_names08[i]
  argument_file2[i,2] <- 2008
  argument_file2[i,3] <- 2008
  argument_file2[i,4] <- 'all'
  argument_file2[i,5] <- 'TRUE'
  argument_file2[i,6] <- 'TRUE'
  argument_file2[i,7] <- 'FALSE'
  argument_file2[i,8] <- 'FALSE'
  argument_file2[i,9] <- 2008
  argument_file2[i,10] <- 'TRUE'
  
  write.csv2(argument_file2[i, ], paste0('BMP_', euring_codes08[i],'_1_63_arg_input_stratum.csv'), row.names = FALSE)
}



##############################################################
#### PREPARING THE TOV-E DATA FOR THE PECBMS REPORT 2023  ####
#### This chunck of code is from Retrieve_Data_from_DB_2023 ##
##############################################################

#folder <- "C:/Users/diego.pavon-jordan/OneDrive - NINA/Documents/PECBMS/2023_Analysis"  
#setwd(folder)


## Packages 
suppressMessages(require('RPostgres'))
suppressMessages(require('odbc'))
suppressMessages(require(RPostgreSQL))
suppressMessages(require(DBI))
suppressMessages(require(tidyverse))
suppressMessages(require(rpostgis))
suppressMessages(require(sf))
suppressMessages(require(lubridate))
suppressMessages(require(xtable))
require('rtrim')


## Connecting to the DB
sort(unique(odbcListDrivers()[[1]]))

con <- DBI::dbConnect(odbc(),
                      Driver   = "SQL server", 
                      Server   = "ninsql07.nina.no",
                      Database = "TOVTaksering",
                      Trusted_Connection = "True")

# Getting entire dataset (incl. -1s and 0s) for all species
Bird_data <- tbl(con, 'TrimResults') %>%
  filter(Year > 2005,
         Species > 0) %>%
  arrange(Site, Species, Year) %>%
  as_tibble() %>%
  select(Id, TrimId, Site,  Year, Species, Count)

unique(Bird_data$Year)

# How many species?
Bird_data %>% select(Species) %>% distinct()
# 230 species

# How many sites?
Bird_data %>% select(Site) %>% distinct()
# 493 sites 


## Linking the Species code from Bird_Data
## with the ArtsID from the Art table in TOVTaksering to get EURING codes

## Getting the data from Art table
SppID <- tbl(con, 'Art') %>%
  mutate(Species = as.numeric(ArtsID)) %>%
  select(Species, EURINGCode, Artsnavn_Lat, FK_Kode_Flokk) %>%
  as_tibble() %>%
  filter(Species > 0) %>%
  arrange(Species) 
# 342 species


## Adding the EURING code form Art to the Bird_data (TrimResults)
Bird_data_EURING <- Bird_data %>% left_join(SppID)

## Do we have missing EURING codes?
Bird_data_EURING %>% filter(is.na(EURINGCode)) %>% distinct(Artsnavn_Lat)
# Yes! 4 species (2979 Columba livia, 9220 Sterna sp., 9280 Loxia sp., 9998)


## I remove these 4 species to have clean data from TrimResults with the EURING code added.
Bird_data_EURING_clean <- Bird_data_EURING %>% filter(!is.na(EURINGCode)) %>%
  mutate(Species_nr = as.numeric(EURINGCode))

## Do we have missing EURING codes?
Bird_data_EURING_clean %>% filter(is.na(EURINGCode))
# No! Good!




## Now, I need to select only the 71 species (as per 2022 list) that we have to report to PECBMS


## Get the selected 71 species and years to be included in the report
## from the PECBMS_2022.R

Species_info <- read.csv2('ReportSpeciesPECBMS2020.csv')

Spp_selection <- read.csv2('PECBMS_species_list_2022.csv') %>%
  left_join(Species_info) %>%
  select(EURINGCode, Year_First, TCF_Basetime, TCF_Slope_From) %>%
  as_tibble()%>%
  mutate(Stratum = 63)
# 71 species - Ok!


Spp_selection <- Spp_selection %>% mutate(Species_nr = EURINGCode)
# 71 species


# Selecting the species from the whole dataset by the EURING code (renamed as Species_nr)
Data_PECBMS <- Bird_data_EURING_clean %>%
  left_join(Spp_selection, by = 'Species_nr')

Data_PECBMS_clean <- Data_PECBMS %>% 
  filter(!is.na(EURINGCode.y))
# 595,051 rows

## Do I have 71 species?
Data_PECBMS_clean %>% distinct(Species) 
# Yes!

write.csv2(Data_PECBMS_clean, 'Selected_species_to_PECBMS_2023.csv', row.names = F)
###############################################################################
## Now in Data_PECBMS_clean we have the full year-site counts             #####
## including -1s and 0s for the 71 species that we report to PECBMS       ##### 
###############################################################################


#### In case you start from here, just load the file :)

## Missing years are already added in the database.
for_trim <- read.csv2('Selected_species_to_PECBMS_2023.csv') %>% 
  select(Site, Year, Species_nr, Count, Artsnavn_Lat) %>%
  as_tibble() %>%
  mutate(EURINGCode = Species_nr)

for_trim %>% select(Species_nr) %>% distinct()
# 71 species. Ok!


## removing the annoying blank space at the end of the names
spp_names <- for_trim %>% 
  select(Artsnavn_Lat) %>% 
  distinct() %>% 
  mutate(Spp_name = str_trim(Artsnavn_Lat))



## make a table indicating the species that need 2007 removed (spp starting in 2008)
euring_codes08_2 <- as.data.frame(euring_codes08) %>%
  mutate(remove_year = 2007) %>%
  as_tibble() %>%
  rename(Species_nr = euring_codes08)

for_trim <- for_trim %>% left_join(spp_names)

## EURING codes for all 71 species
euring_codes <- for_trim %>%
  select(EURINGCode) %>%
  distinct() %>%
  arrange(EURINGCode)


# Create a list for all species.
# The new trim does not like -1s, so we need to add NAs instead
# I also remove here the undesired years
by_spp <- for_trim %>%
  rename(year = Year,
         site = Site,
         count = Count) %>%
  select(-Artsnavn_Lat) %>%
  #na_if(-1) %>%
  filter(year > 2006) %>%
  left_join(euring_codes08_2) %>%
  arrange(site, year) %>%
  mutate(remove_row = case_when(year == remove_year ~ 'yes',
                                TRUE ~'no')) %>%
  filter(remove_row == 'no') %>%
  select(site, year, count, EURINGCode) %>%
  arrange(EURINGCode) %>%
  group_by(EURINGCode) %>%
  nest()

by_spp

#check it out
by_spp[['data']][[1]] %>% filter(site == 101) # starts in 2008!
euring_codes[1,]
filter(for_trim, EURINGCode == 1860)
# OK!


## Saving the species-specific input files as per rtrim-shell format request
for(i in 1:nrow(euring_codes)){
  write.csv2(by_spp[['data']][[i]], paste0('BMP_', euring_codes[i,1],'_1_63_counts.csv'), row.names = FALSE)
}


################################################################################
################# END OF DATA PREPARATION ######################################
################################################################################




################################################################################
#### START CODE PROVIDED BY PECBMS #############################################
#### Running rtrim_shell           #############################################
#### Version RTRIM-shell_v1.4      #############################################
################################################################################

################################################################################
# PREPARATION.
################################################################################
# 1-06-2020 Eva Silarova unified upper and lower cases in variable names - v1.1
#rm(list=ls())

# Link to script with functions   !! user dependent !!
source("C:/Users/diego.pavon-jordan/Documents/PECBMS/02_RTRIM_shell_v1.4/functions_Shell_Rtrim.R")

# The required library
#library(rtrim)

# Select the directory which contains the RTRIM output files
# This directory also contains the files with arguments used in the rtrim function.
# folder <- choose.dir(default ="D:/WORK/RTRIM-shell/Examples_for_coordinators_Input-files/00_RTRIM-shell_inputs&scripts_corrected_AG", caption = "Select directory which contains data which has been analysed")
setwd(folder)

# Select all files in the directory which contain arguments.
# The files can be recognized by the pattern "arg_input_arguments".
listSpeciesStratumCombinations <- dir(folder, pattern = "arg_input_stratum.csv")

# The number of files in the directory.                                                
numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# Make a table for an overview. 
# The tabel lists which analyses were successful or not for all runs 

overview <- makeOverview(listSpeciesStratumCombinations)
# The code of this function can be found in functions_Shell_Rtrim.R.
# The source command (see above) makes it possible to use this function in this code.

# --> ALL WORKING GOOD!




#####################################################################################################################
# READING ARGUMENT AND COUNT FILES. 
#####################################################################################################################

for (j in 1:numberSpeciesStratumCombinations) {
  
  # j<- 1
  
  # The file with arguments contains the information to analyse the counts for a particular combination of species and stratum. 
  # The arguments are used when calling the function 'rtrim'
  arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
  
  # The file with counts contains the counts for a particular combination of species and stratum.
  # Weights may also be present in this file.
  counts <- read.csv2(paste(arguments$File, "_counts.csv", sep = ""), stringsAsFactors = FALSE, header = TRUE)
  
  overview$first_year[j] <- min(counts$year, na.rm = TRUE)
  overview$last_year[j]  <- max(counts$year, na.rm = TRUE)
  # First and last year for each combination of stratum and species is stored for later use.
  
  #####################################################################################################################
  # RUNNING RTRIM.
  #####################################################################################################################
  
  # Start with the most elaborate model and switch automatically to a more simple model when needed.
  
  result <- tryCatch(
    {
      # Attempt 1
      if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
        trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
      } else {
        if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
          trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        } else{
          if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
            trim(count ~ site + (year + month), data = counts,            model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
          } else{
            trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
          }
        }
      }
    }
    , error = warning)     
  
  if (class(result) == "trim") {
    
    save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
    overview$attempt_1[overview$ss_combinations == arguments$File] <- "success"
    overview$success[overview$ss_combinations == arguments$File] <- "yes"
    
  } else {
    
    # First attempt failed, try a less elaborate model by setting serial correlation off. 
    # Also, when month factors are available, no changepoints are estimated in the next model.
    
    overview$attempt_1[overview$ss_combinations == arguments$File] <- "error"
    overview$error_1[overview$ss_combinations == arguments$File] <- result
    
    result <- tryCatch(
      {
        # attempt 2 
        if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
          trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        } else {
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
            trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          } else{
            if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
              trim(count ~ site + (year + month), data = counts,                      model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
            } else{
              trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
            }
          }
        }
      }
      , error = warning)
    
    if (class(result) == "trim") {
      
      save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
      overview$attempt_2[overview$ss_combinations == arguments$File] <- "success"
      overview$success[overview$ss_combinations == arguments$File] <- "yes"
      
    } else {
      
      # Second attempt also failed. Now try an even more simple model: no month factors, no changepoints, but serial correlation switched on. 
      # Note that no further options are available to include month factors in the model. 
      
      overview$attempt_2[overview$ss_combinations == arguments$File] <- "error"
      overview$error_2[overview$ss_combinations == arguments$File] <- result
      
      if (arguments$Presence_monthfactors == TRUE) {
        
        cat("Analysis failed for this combination of species and stratum:", arguments$File, "\n")
        
      }
      
      result <- tryCatch( 
        {
          # attempt 3
          
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
            
            trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
            
          } else{
            
            if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
              trim(count ~ site + year, data = counts,                      model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
            }
          }
        }
        , error = warning)
      
      if (class(result) == "trim") {
        
        save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
        overview$attempt_3[overview$ss_combinations == arguments$File] <- "success"
        overview$success[overview$ss_combinations == arguments$File] <- "yes"
        
      } else {
        if (arguments$Presence_monthfactors == FALSE) {
          overview$attempt_3[overview$ss_combinations == arguments$File] <- "error"
          overview$error_3[overview$ss_combinations == arguments$File] <- result
        }
        
        # Third attempt also failed. Now try the most simple model: no changepoints at all and no serial correlation. 
        
        result <- tryCatch(
          {
            # Final attempt
            
            if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
              trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
            } else{
              if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                trim(count ~ site + year, data = counts,                      model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              } 
            }
            
          }
          
          , error = warning)
        
        if (class(result) == "trim") {
          
          save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
          overview$attempt_4[overview$ss_combinations == arguments$File] <- "success"
          overview$success[overview$ss_combinations == arguments$File] <- "yes"
          
        } else {
          if (arguments$Presence_monthfactors == FALSE) {
            
            # When this analysis also fails, send a error message to screen.
            
            overview$attempt_4[overview$ss_combinations == arguments$File] <- "error"
            overview$error_4[overview$ss_combinations == arguments$File] <- result
            
            cat("Analysis failed for this combination of species and stratum:", arguments$File, "\n")
          }        
        }
      }  
    } 
  }  
}

#####################################################################################################################
# WRITING OVERVIEW OF RTRIM SUCCESSES and FAILURES.  
#####################################################################################################################
overview <- overview[order(overview$species_number, overview$stratum_number), ]
write.csv2(overview, "overview.csv", row.names = FALSE)





#####################################################################################################################
# RTRIM-shell second script. Version RTRIM-shell_v1.3
# Tool to process rtrim results (= TRIM in R) for multiple species and subsets of sites.
# Marnix de Zeeuw rtrim@cbs.nl Statistics Netherlands 2019.
#####################################################################################################################

#####################################################################################################################
# SELECTING SUCCESSFUL RUNS.
#####################################################################################################################
# Determine which datafiles (combinations of species and stratum) have been analysed successfully. 
# This information is found in "overview.csv", which has been produced by the script that called rtrim.

overview <- read.csv2("overview.csv")
listsuccessfulAnalyses <- overview$ss_combinations[overview$success == "yes"]
listSpeciesStratumCombinations <- paste(listsuccessfulAnalyses, "_arg_input_stratum.csv", sep = "")

# Determine how many combinations of species and stratum have been analysed successfully.
numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# File with trends and indices for each combination of species and stratum
all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)

#####################################################################################################################
# PROCESSING OUTPUT OF SUCCESFUL RUNS.
#####################################################################################################################

for (j in 1:numberSpeciesStratumCombinations) {
  
  # The file with arguments contains the arguments to run the analysis for a particular combination of species and stratum (stratum is e.g. a region). 
  # These arguments are used when calling the function 'rtrim'.
  arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
  counts <- read.csv2(paste(arguments$File, "_counts.csv", sep = ""), header = TRUE)
  
  # Also the result of the rtrim-analysis is required.
  load(paste(arguments$File, ".RData", sep = "")) # produces object with name "result"
  
  #####################################################################################################################
  # CREATING FILES FOR LATER USE. 
  #####################################################################################################################
  
  # Several output files are created (dataframes). Filling the files with output is done at a later stage.
  # To create the file containing indices and time totals (indices_TT_file) requires the file with results.
  # The file arg_output contains the slopes and arguments used to run the rtrim function.
  
  indices_TT_file <- make_Indices_TT_file(result = result)
  
  arg_output_file <- make_arg_output_file(arguments = arguments)
  
  #####################################################################################################################
  # FILLING OUTPUT.
  # A few functions are used to fill the output files.
  # These functions can be found in "functions_Shell_Rtrim.r".
  # This script sources to that file to include them here, which makes it easier to see what is done in this script. 
  #####################################################################################################################
  
  indices_TT_file <- fill_Indices_TT_file(indices_TT_file = indices_TT_file, result = result, arguments = arguments)
  
  arg_output_file <- fill_arg_output_file(arg_output_file = arg_output_file, result = result, counts = counts)
  
  all_Indices_All_Trends <- fill_All_Indices_All_Trends(all_Indices_All_Trends = all_Indices_All_Trends, result = result, arguments = arguments, j = j, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  covariant_matrix <- vcov(result)
  
  if (arguments$Save_fitted_values){
    
    FI <- results(result)
    
  }
  
  #####################################################################################################################
  # WRITING OUTPUT FILES. 
  #####################################################################################################################
  # Indices and time totals.
  name_Indices_TT_file <- paste(arguments$File, "_indices_TT.csv", sep = "")
  write.csv2(indices_TT_file, name_Indices_TT_file, row.names = FALSE)
  ############################################################ 
  # Slopes entire period and arguments.
  name_arg_output_file <- paste(arguments$File, "_arg_output.csv", sep = "") 
  write.csv2(arg_output_file, name_arg_output_file, row.names = FALSE)
  ############################################################
  # Covariant matrix.
  name_covariant_matrix <- paste(arguments$File, "_ocv.csv", sep = "")
  write.csv2(covariant_matrix, name_covariant_matrix, row.names = FALSE)
  ############################################################
  # File with fitted values.
  if (arguments$Save_fitted_values){
    name_Fitted_Values_File <- paste(arguments$File, "_fitted_values.csv", sep = "")
    write.csv2(FI, name_Fitted_Values_File, row.names = FALSE)
  }
  ############################################################
  #25.01.2022 Javi: Generating Tables 1-3 for the Coordinators
  name_rawdata=paste(arguments$File, "_summarizing_tables.txt", sep = "")
  #We are filling the  table on numbers of sites and counts
  write(paste0("\n", arguments$File, "\n"), file = name_rawdata)
  
  Table1=data.frame(Parameter=c("Sites","Times","Number of observerd zero counts","Number of observed positive counts",
                                "Total number of observed counts", "Number of missing counts", "Total number of counts"),
                    Value= c(arg_output_file$N_sites,arg_output_file$N_time_values,arg_output_file$N_observed_zero_counts,
                             arg_output_file$N_observed_positive_counts,(arg_output_file$N_observed_positive_counts+arg_output_file$N_observed_zero_counts),
                             arg_output_file$N_missing_counts,arg_output_file$N_counts_total))
  
  
  
  #Table showing those sites containing more than the 10% of observations
  
  ?data.frame
  total_count=sum(counts$count[!is.na(counts$count)])
  tenpercent=(total_count*10)/100
  
  Table2=as.data.frame(matrix(NA,ncol = 3))
  names(Table2)=c("Site Number","Observed total","%")
  xx=subset(counts,!is.na(counts$count)==T)
  y=unique(xx$site)
  for(i in 1:length(unique(xx$site))){
    x=subset(xx,xx$site==y[i])
    s=sum(x$count)
    if (s>tenpercent){
      if(is.na(Table2[1,1])==T){
        Table2$`Site Number`=x$site[1]
        Table2$`Observed total`=s
        percentage=(s*100/total_count)
        Table2$`%`=percentage
      }else{
        percentage=(s*100/total_count)
        vector=c(x$site[1],s,percentage)
        Table2=rbind(Table2,vector)
      }
    }
  }  
  
  ##Table 3 is timpoint average table
  
  Table3=as.data.frame(matrix(NA,ncol=5,nrow=length(unique(counts$year[!is.na(counts$count)]))))
  
  names(Table3)=c("Timepoint","Observations","Average","Index","Real_Number")
  years=unique(counts$year)
  for(i in 1:length(years)){
    Table3$Timepoint[i]=years[i]
    visitedsites=subset(xx,xx$year==years[i])
    Table3$Observations[i]=nrow(visitedsites)
    Table3$Real_Number[i]=sum(visitedsites$count)
    Table3$Average[i]=Table3$Real_Number[i]/Table3$Observations[i]
    if(i==1){
      Table3$Index[i]=1.00
    }else{
      Table3$Index[i]=Table3$Average[i]/Table3$Average[1]
    }
  } 
  #This is the last calculation for Table 1 that requires data from table 3
  x=0
  for(i in 1:nrow(Table3)){
    s=Table3$Observations[i]*Table3$Average[i]
    x=x+s
  }
  
  Total_count=c("Total count",x)
  Total_count
  Table1=rbind(Table1,Total_count)
  #Filling the txt!
  write("\n 1.	Summarizing table on numbers of sites and counts \n", file = name_rawdata, append = T)
  capture.output( print.data.frame(Table1, row.names=F, print.gap=3, quote=F, right=F),  file=name_rawdata,append = T  )
  write("\n 2.	Sites containing more than 10% of the total count \n", file = name_rawdata, append = T)
  capture.output( print.data.frame(Table2, row.names=F, print.gap=3, quote=F, right=F),  file=name_rawdata,append = T  )
  write("\n 3.	Time Point Averages table \n", file = name_rawdata, append = T)
  capture.output( print.data.frame(Table3, row.names=F, print.gap=3, quote=F, right=F),  file=name_rawdata,append = T  )
}

all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
write.csv2(all_Indices_All_Trends, "All_Indices_All_Trends.csv", row.names = FALSE)           # 09-04-2021 Eva Silarova all changed to capital All: "all_Indices_All_Trends.csv" was changed to "All_Indices_All_Trends.csv" - v1.3
############################################################




## Changing the name of the files (removing the 'BMP_' part)
## first copy the folder with all results and use the new folder. Set it as wd
## Making the new folder is optional. You can just use the same wd

# create new folder
new_folder <- 'Species_files'
dir.create(file.path(folder, new_folder))

# copy everything in that folder
species_files <- dir(folder, pattern = ".csv")
sumar_tables <- dir(folder, pattern = ".txt")

# custom function
my_function <- function(x){
  file.copy(from = file.path("C:/Users/diego.pavon-jordan/OneDrive - NINA/Documents/PECBMS/2023_Analysis", x) ,
            to = file.path("C:/Users/diego.pavon-jordan/OneDrive - NINA/Documents/PECBMS/2023_Analysis/Species_files", x) )
}

# apply the function to all files
lapply(species_files, my_function)
lapply(sumar_tables, my_function)


setwd(paste0(folder, '/Species_files'))

Old_names <- list.files()

New_names <- list.files() %>%
  str_replace('BMP_', '')

file.rename(Old_names, New_names)





















