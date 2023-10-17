# RSWAN Version 6 (script 03_Swan_EuromonitoringV6.R=v5, but Pools and Schedules are changed according to changes in the database)
# Arco van Strien - Okt 2020
# Last corrections were done 16-2-2021

# New in version 2
# Select all available data files of Counties (Level1 files) (using new function Select_all_available_datafiles). This is optional: if B_sel=0 in function list, this is skipped 
# Swan_schedule automatically sorted to ensure proper data processing
# Output folder added for result files
# Output files with all results in vertical format added, both for countries and pools
# Output files with al results in horizontal format changed, both for countries and pools
# Output files of results of countries and pools combined
# Output files of individual species and stratum converted into format similar as output from RTRIM shell (using new function Make_empty_arg_output_file)
# When base year does not occur in data, first available year in data treated as base year
# Function truncation of years added, using year_first and year_last as specified in Species-Country file
# Error fixed in range of years (using new function Assess_range_of_years_for_resultsfile)


# Start of code: 


Check_steering_files <- function(general_folder_path) {
  ###########################################################################
  # Aim: Check presence and consistency of steering files                   #
  # Test if all steering and lookup files occur in the general folder       # 
  # Run this function before any other program                              #
  # Fill in the path of the general folder in the call of this function     #
  ###########################################################################
  
  setwd(general_folder_path) 
  
  cat("\n")
  cat("@@@ function Check_steering_files executed", "\n")
  cat("\n")
  
  ###########################################################################
  # Check existence and consistency individual steering files               #                   
  ###########################################################################
  
  if (file.exists("Schedules.csv")  == FALSE)    { cat("Error: missing file Schedules.csv","\n")    }  
  if (file.exists("Schedules.csv")  == TRUE)     { 
    #   check for double schedule codes
    schedules_file                  <- read.table("Schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    check_double_schedules          <- schedules_file[, c("Schedule_code")]  
    a <- length(check_double_schedules)
    b <- length(unique(check_double_schedules)) 
    if (a != b) { cat("Error: the same schedule code occurs multiple times in Schedules.csv", "\n") }
  } # end of check Schedules.csv
  
  
  if (file.exists("Pools.csv")      == FALSE)    { cat("Error: missing file Pools.csv","\n")        }  
  if (file.exists("Pools.csv")      == TRUE)     { 
    #   check for double pool codes
    pools_file                      <- read.table("Pools.csv",header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    check_double_pools              <- pools_file[, c("Pool")]  
    a <- length(check_double_pools)
    b <- length(unique(check_double_pools)) 
    if (a != b) { cat("Error: the same pool occurs multiple times in Pools.csv", "\n")              }
  } # end of check Pools.csv
  
  
  if (file.exists("Countries.csv")  == FALSE)    { cat("Error: missing file Countries.csv","\n")    }  
  if (file.exists("Countries.csv")  == TRUE)     { 
    #   check for double country codes
    countries_file                  <- read.table("Countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    check_double_countries          <- countries_file[, c("Country_code")] 
    a <- length(check_double_countries)
    b <- length(unique(check_double_countries)) 
    if (a != b) { cat("Error: the same country occurs multiple times in Countries.csv", "\n")       }
  } # end of check Countries.csv
  
  
  if (file.exists("Species.csv")    == FALSE)    { cat("Error: missing file Species.csv","\n")      }  
  if (file.exists("Species.csv")    == TRUE)     { 
    #   check for double species codes 
    species_file                    <- read.table("Species.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol     
    check_double_species            <- species_file[, c("Species_nr")]  
    a <- length(check_double_species)
    b <- length(unique(check_double_species)) 
    if (a != b) { cat("Error: the same species occurs multiple times in Species.csv", "\n")         }
  } # end of check Species.csv
  
  
  if (file.exists("Species_countries.csv")  == FALSE) { cat("Error: missing file Species_countries.csv","\n")    }  
  if (file.exists("Species_countries.csv")  == TRUE)  { 
    #   check for double species_country codes
    species_countries_file                  <- read.table("Species_countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol
    check_double_species_countries          <- species_countries_file[, c("Species_nr", "Country_code", "Population size (geomean)")]  # population size is needed below
    check_double_species_countries$combination_column <-paste(as.character(check_double_species_countries[,"Species_nr"]), as.character(check_double_species_countries[,"Country_code"]), sep="-") 
    a <- nrow(check_double_species_countries)
    b <- length(unique(check_double_species_countries[, "combination_column"])) 
    if (a != b) { cat("Error: the same country occurs multiple times in Species_countries.csv", "\n") }
    if (!is.numeric(species_countries_file$`Population size (geomean)`)==T){ #Javi 16/02/2023 Aditional controls to check if formating is correct
      cat("Error: Species' population size is not numeric in Species_countries.csv", "\n")}
    if (!is.numeric(species_countries_file$`Population size_max`)==T){ #Javi 16/02/2023 Aditional controls to check if formating is correct
       cat("Error: Species' population size max is not numeric in Species_countries.csv", "\n")
    }# End of Javi 16/02/2023
    
  } # end of check Species_countries.csv
  
  ###########################################################################
  # Check mutual consistency of steering files                              #                   
  ###########################################################################
  
  if (file.exists("Swan_schedules.csv")  == FALSE)    { cat("Error: missing file Swan_schedules.csv","\n")   }  
  if (file.exists("Swan_schedules.csv")  == TRUE)     { 
   
      Swan_schedules  <- read.table("Swan_schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    
    
    #   test if all schedules exist in file with schedules
    schedules_in_Swan_schedules <- as.character(unique(Swan_schedules[, "Schedule_code"]))
    schedules_found <-  is.element(schedules_in_Swan_schedules, schedules_file[, "Schedule_code"])  
    if (length(schedules_in_Swan_schedules) != sum(schedules_found)) { # some schedules are lacking
      missing_schedule <- which(schedules_found == FALSE) 
      cat("Error: Schedule", schedules_in_Swan_schedules[missing_schedule], "does not occur in Schedules.csv", "\n") 
    }
    
    
    #   test if all pools exist in file with schedules
    LevelCols <- as.vector(grep("Level", names(Swan_schedules)))
    for (k in 2:length(LevelCols)) { # first level column contains no pools but countries
      pools_in_Swan_schedules <- unique(Swan_schedules[, LevelCols[k]])
      pools_found <-  is.element(pools_in_Swan_schedules, pools_file[, "Pool"])  
      if (length(pools_in_Swan_schedules) != sum(pools_found)) { # some pools are lacking
        missing_pool <- which(pools_found == FALSE) 
        cat("Error: Pool", pools_in_Swan_schedules[missing_pool], "in level", k, "of Swan_schedules.csv does not occur in Pools.csv", "\n") 
      }
    }
    
    
    #   test if all countries exist in file with schedules
    countries_in_Swan_schedules <- unique(Swan_schedules[, "Level1"])
    countries_found <-  is.element(countries_in_Swan_schedules, countries_file[, "Country_code"])  
    if (length(countries_in_Swan_schedules) != sum(countries_found)) { # some countries are lacking
      missing_country <- which(countries_found == FALSE) 
      cat("Error: Country", countries_in_Swan_schedules[missing_country], "does not occur in Countries.csv", "\n") 
    }
    
    
    #   test if all species exist in file with schedules
    species_in_Swan_schedules <- unique(Swan_schedules[, "Species_nr"])
    species_found <-  is.element(species_in_Swan_schedules, species_file[, "Species_nr"])  
    if (length(species_in_Swan_schedules) != sum(species_found)) { # some species are lacking
      missing_species <- which(species_found == FALSE) 
      cat("Error: Species", species_in_Swan_schedules[missing_species], "does not occur in Species.csv", "\n") 
    }
    
    
    #   test if weights deviate from 1 (if so, no problem)
    nr_weight_1 <- sum(Swan_schedules[, "Weight"] == 1)
    if (nr_weight_1 == nrow(Swan_schedules)) { cat("All weights in Swan_schedule.csv have default value 1", "\n") }
    if (nr_weight_1 < nrow(Swan_schedules))  { cat("Note: some weights in Swan_schedule.csv have value != default value  1", "\n") }
    
    ###########################################################################
    #   Check consistency of Swan_schedules.csv                               #                   
    ###########################################################################
    
    #   a species_country combination cannot be attributed to more than 1 pool
    Swan_schedules$CombiCol1 <- paste(as.character(Swan_schedules$Schedule), as.character(Swan_schedules$Species), as.character(Swan_schedules$Level1),sep="-") 
    no_of_unique_records <- length(unique(Swan_schedules$CombiCol1))
    if (no_of_unique_records != nrow(Swan_schedules)) { cat("Error: one or more countries have been attributed to multiple pools in level2 in Swan_schedules.csv", "\n") }
    
    
    #   the same species_pool combination cannot be attributed to more than 1 pool of a higher level
    LevelCols <- as.vector(grep("Level", names(Swan_schedules))) # note: poolno's are allowed to overlap with country no's
    for (k in 2:length(LevelCols)) { 
      Swan_schedules$CombiCol2 <- paste(as.character(Swan_schedules$Schedule), as.character(Swan_schedules$Species), as.character(Swan_schedules[, LevelCols[k-1]]), sep="_") 
      Swan_schedules$CombiCol3 <- paste(as.character(Swan_schedules$Schedule), as.character(Swan_schedules$Species), as.character(Swan_schedules[, LevelCols[k]]), sep="_") 
      #     inspect for each level whether CombiCol2 conflicts with CombiCol3
      pools_per_level <- unique(Swan_schedules[, "CombiCol2"])
      for (p in 1: length(pools_per_level)) { 
        pools_in_combi3 <- unique(Swan_schedules[Swan_schedules$CombiCol2 == pools_per_level[p], "CombiCol3"])
        if (length(pools_in_combi3) > 1) { cat("Error: Pool", pools_per_level[p], "in lower level used in multiple pools of higher level:", pools_in_combi3, "\n") }
      }      
    }
    
    
    #   test existence Swan_function_list and also if there are functions selected in each selected schedule
    if (file.exists("Swan_function_list.csv")== FALSE) { cat("Error: missing file Swan_function_list.csv","\n")   }  
    if (file.exists("Swan_function_list.csv")== TRUE)  { 
      Swan_function_list     <- read.table("Swan_function_list.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")               # function list
      schedules_selected_in_function_list  <- as.character(unique(Swan_function_list[Swan_function_list$B_sel == 1,"Schedule_code"])) # B_sel = 1 means selected
      schedules_selected_in_Swan_schedules <- as.character(unique(Swan_schedules[Swan_schedules$B_sel == 1,"Schedule_code"]))     
      schedules_selected_found <-  is.element(schedules_selected_in_Swan_schedules, schedules_selected_in_function_list)  
      if (sum(schedules_selected_found) !=  length(schedules_selected_in_Swan_schedules) ) { 
        cat("Error: some schedules selected in Swan_schedules.csv have no functions selected in Swan_function_list.csv", "\n") 
      }  
      
    } # end if Swan_function_list exists
    
    #   test if any species*countries are selected in Swan_schedules
    species_countries_selected <- Swan_schedules[, "B_sel"] == 1
    if (sum(species_countries_selected) == 0) { cat("Error! No species*countries selected in Swan_schedules.csv", "\n")  }
    
    #   ensure proper sorting of Swan_schedules.csv = sorting from highest Level to Lowest level  
    nlevel <- length(LevelCols)  # take into account that no of Levels might differ from the default 4 levels 
    
    if (nlevel == 3)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if (nlevel == 4)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"],Swan_schedules[,"Level4"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if (nlevel == 5)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"],Swan_schedules[,"Level5"], Swan_schedules[,"Level4"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if (nlevel == 6)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"],Swan_schedules[,"Level6"],Swan_schedules[,"Level5"], Swan_schedules[,"Level4"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if ((nlevel < 3) | (nlevel > 6)) { cat("Warning: ensure that the sorting of Swan_schedules.csv is adequate", "\n") }
    
    Swan_schedules$CombiCol1	<- NULL
    Swan_schedules$CombiCol2	<- NULL
    Swan_schedules$CombiCol3	<- NULL
    write.table(Swan_schedules, file = "Swan_schedules.csv",  dec=".",sep=";",row.names=FALSE) # properly sorted #13/02/2023 JAVI, DECIMAL SYMBOL
    
  } # end if Swan_schedules exists
  
  #Check for missing species-country combinations in species_countries.csv that are present in the schedule.csv #Javi 16/02/2023 new controls to prevent errors in rswan execution
  
  swan_schedule=read.table("Swan_schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
  Swan_schedules_selected <- Swan_schedules[Swan_schedules$B_sel == 1,]  
  speciescountry_code<-paste(species_countries_file$Species_nr,species_countries_file$Country_code)
  schedule_speciescountry <- paste(Swan_schedules_selected$Species_nr,Swan_schedules_selected$Level1)
  
  missing_combinations<-schedule_speciescountry[!(schedule_speciescountry %in% speciescountry_code) ]
  if(length(missing_combinations)!=0){
    missing_combinations <- gsub(" ","-",missing_combinations)
    cat("Error: The following species-country combination are selected in the schedule but not present in the species country file","\n" )
    cat("\n")
    cat(missing_combinations,"\n")
  } # END Check for missing species-country combinations # Javi 16/02/2023
  
  
} # E N D of function Check_steering_files
##########################################################################################################################################################
#20/04/2023 Javi: I duplicated the function check steering files and adapted to perform their checks to the combine files
Check_steering_files_combine <- function(general_folder_path) {
  ###########################################################################
  # Aim: Check presence and consistency of steering files                   #
  # Test if all steering and lookup files occur in the general folder       # 
  # Run this function before any other program                              #
  # Fill in the path of the general folder in the call of this function     #
  ###########################################################################
  
  setwd(general_folder_path) 
  
  cat("\n")
  cat("@@@ function Check_steering_files executed", "\n")
  cat("\n")
  
  ###########################################################################
  # Check existence and consistency individual steering files               #                   
  ###########################################################################
  
  if (file.exists("Schedules.csv")  == FALSE)    { cat("Error: missing file Schedules.csv","\n")    }  
  if (file.exists("Schedules.csv")  == TRUE)     { 
    #   check for double schedule codes
    schedules_file                  <- read.table("Schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    check_double_schedules          <- schedules_file[, c("Schedule_code")]  
    a <- length(check_double_schedules)
    b <- length(unique(check_double_schedules)) 
    if (a != b) { cat("Error: the same schedule code occurs multiple times in Schedules.csv", "\n") }
  } # end of check Schedules.csv
  
  
  if (file.exists("Pools.csv")      == FALSE)    { cat("Error: missing file Pools.csv","\n")        }  
  if (file.exists("Pools.csv")      == TRUE)     { 
    #   check for double pool codes
    pools_file                      <- read.table("Pools.csv",header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    check_double_pools              <- pools_file[, c("Pool")]  
    a <- length(check_double_pools)
    b <- length(unique(check_double_pools)) 
    if (a != b) { cat("Error: the same pool occurs multiple times in Pools.csv", "\n")              }
  } # end of check Pools.csv
  
  
  if (file.exists("Countries.csv")  == FALSE)    { cat("Error: missing file Countries.csv","\n")    }  
  if (file.exists("Countries.csv")  == TRUE)     { 
    #   check for double country codes
    countries_file                  <- read.table("Countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    check_double_countries          <- countries_file[, c("Country_code")] 
    a <- length(check_double_countries)
    b <- length(unique(check_double_countries)) 
    if (a != b) { cat("Error: the same country occurs multiple times in Countries.csv", "\n")       }
  } # end of check Countries.csv
  
  
  if (file.exists("Species.csv")    == FALSE)    { cat("Error: missing file Species.csv","\n")      }  
  if (file.exists("Species.csv")    == TRUE)     { 
    #   check for double species codes 
    species_file                    <- read.table("Species.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol     
    check_double_species            <- species_file[, c("Species_nr")]  
    a <- length(check_double_species)
    b <- length(unique(check_double_species)) 
    if (a != b) { cat("Error: the same species occurs multiple times in Species.csv", "\n")         }
  } # end of check Species.csv
  
  
  if (file.exists("Species_countries.csv")  == FALSE) { cat("Error: missing file Species_countries.csv","\n")    }  
  if (file.exists("Species_countries.csv")  == TRUE)  { 
    #   check for double species_country codes
    species_countries_file                  <- read.table("Species_countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol
    check_double_species_countries          <- species_countries_file[, c("Species_nr", "Country_code", "Population size (geomean)")]  # population size is needed below
    check_double_species_countries$combination_column <-paste(as.character(check_double_species_countries[,"Species_nr"]), as.character(check_double_species_countries[,"Country_code"]), sep="-") 
    a <- nrow(check_double_species_countries)
    b <- length(unique(check_double_species_countries[, "combination_column"])) 
    if (a != b) { cat("Error: the same country occurs multiple times in Species_countries.csv", "\n") }
    if (!is.numeric(species_countries_file$`Population size (geomean)`)==T){ #Javi 16/02/2023 Aditional controls to check if formating is correct
      cat("Error: Species' population size is not numeric in Species_countries.csv", "\n")}
    if (!is.numeric(species_countries_file$`Population size_max`)==T){ #Javi 16/02/2023 Aditional controls to check if formating is correct
      cat("Error: Species' population size max is not numeric in Species_countries.csv", "\n")
    }# End of Javi 16/02/2023
    
  } # end of check Species_countries.csv
  
  ###########################################################################
  # Check mutual consistency of steering files                              #                   
  ###########################################################################
  
  if (file.exists("Swan_schedules_combine.csv")  == FALSE)    { cat("Error: missing file Swan_schedules.csv","\n")   }  
  if (file.exists("Swan_schedules_combine.csv")  == TRUE)     { 
    
    Swan_schedules  <- read.table("Swan_schedules_combine.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    
    
    #   test if all schedules exist in file with schedules
    schedules_in_Swan_schedules <- as.character(unique(Swan_schedules[, "Schedule_code"]))
    schedules_found <-  is.element(schedules_in_Swan_schedules, schedules_file[, "Schedule_code"])  
    if (length(schedules_in_Swan_schedules) != sum(schedules_found)) { # some schedules are lacking
      missing_schedule <- which(schedules_found == FALSE) 
      cat("Error: Schedule", schedules_in_Swan_schedules[missing_schedule], "does not occur in Schedules.csv", "\n") 
    }
    
    
    #   test if all pools exist in file with schedules
    LevelCols <- as.vector(grep("Level", names(Swan_schedules)))
    for (k in 2:length(LevelCols)) { # first level column contains no pools but countries
      pools_in_Swan_schedules <- unique(Swan_schedules[, LevelCols[k]])
      pools_found <-  is.element(pools_in_Swan_schedules, pools_file[, "Pool"])  
      if (length(pools_in_Swan_schedules) != sum(pools_found)) { # some pools are lacking
        missing_pool <- which(pools_found == FALSE) 
        cat("Error: Pool", pools_in_Swan_schedules[missing_pool], "in level", k, "of Swan_schedules.csv does not occur in Pools.csv", "\n") 
      }
    }
    
    
    #   test if all countries exist in file with schedules
    countries_in_Swan_schedules <- unique(Swan_schedules[, "Level1"])
    countries_found <-  is.element(countries_in_Swan_schedules, countries_file[, "Country_code"])  
    if (length(countries_in_Swan_schedules) != sum(countries_found)) { # some countries are lacking
      missing_country <- which(countries_found == FALSE) 
      cat("Error: Country", countries_in_Swan_schedules[missing_country], "does not occur in Countries.csv", "\n") 
    }
    
    
    #   test if all species exist in file with schedules
    species_in_Swan_schedules <- unique(Swan_schedules[, "Species_nr"])
    species_found <-  is.element(species_in_Swan_schedules, species_file[, "Species_nr"])  
    if (length(species_in_Swan_schedules) != sum(species_found)) { # some species are lacking
      missing_species <- which(species_found == FALSE) 
      cat("Error: Species", species_in_Swan_schedules[missing_species], "does not occur in Species.csv", "\n") 
    }
    
    
    #   test if weights deviate from 1 (if so, no problem)
    nr_weight_1 <- sum(Swan_schedules[, "Weight"] == 1)
    if (nr_weight_1 == nrow(Swan_schedules)) { cat("All weights in Swan_schedule.csv have default value 1", "\n") }
    if (nr_weight_1 < nrow(Swan_schedules))  { cat("Note: some weights in Swan_schedule.csv have value != default value  1", "\n") }
    
    ###########################################################################
    #   Check consistency of Swan_schedules.csv                               #                   
    ###########################################################################
    
    #   a species_country combination cannot be attributed to more than 1 pool
    Swan_schedules$CombiCol1 <- paste(as.character(Swan_schedules$Schedule), as.character(Swan_schedules$Species), as.character(Swan_schedules$Level1),sep="-") 
    no_of_unique_records <- length(unique(Swan_schedules$CombiCol1))
    if (no_of_unique_records != nrow(Swan_schedules)) { cat("Error: one or more countries have been attributed to multiple pools in level2 in Swan_schedules.csv", "\n") }
    
    
    #   the same species_pool combination cannot be attributed to more than 1 pool of a higher level
    LevelCols <- as.vector(grep("Level", names(Swan_schedules))) # note: poolno's are allowed to overlap with country no's
    for (k in 2:length(LevelCols)) { 
      Swan_schedules$CombiCol2 <- paste(as.character(Swan_schedules$Schedule), as.character(Swan_schedules$Species), as.character(Swan_schedules[, LevelCols[k-1]]), sep="_") 
      Swan_schedules$CombiCol3 <- paste(as.character(Swan_schedules$Schedule), as.character(Swan_schedules$Species), as.character(Swan_schedules[, LevelCols[k]]), sep="_") 
      #     inspect for each level whether CombiCol2 conflicts with CombiCol3
      pools_per_level <- unique(Swan_schedules[, "CombiCol2"])
      for (p in 1: length(pools_per_level)) { 
        pools_in_combi3 <- unique(Swan_schedules[Swan_schedules$CombiCol2 == pools_per_level[p], "CombiCol3"])
        if (length(pools_in_combi3) > 1) { cat("Error: Pool", pools_per_level[p], "in lower level used in multiple pools of higher level:", pools_in_combi3, "\n") }
      }      
    }
    
    
    #   test existence Swan_function_list and also if there are functions selected in each selected schedule
    if (file.exists("Swan_function_list.csv")== FALSE) { cat("Error: missing file Swan_function_list.csv","\n")   }  
    if (file.exists("Swan_function_list.csv")== TRUE)  { 
      Swan_function_list     <- read.table("Swan_function_list.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")               # function list
      schedules_selected_in_function_list  <- as.character(unique(Swan_function_list[Swan_function_list$B_sel == 1,"Schedule_code"])) # B_sel = 1 means selected
      schedules_selected_in_Swan_schedules <- as.character(unique(Swan_schedules[Swan_schedules$B_sel == 1,"Schedule_code"]))     
      schedules_selected_found <-  is.element(schedules_selected_in_Swan_schedules, schedules_selected_in_function_list)  
      if (sum(schedules_selected_found) !=  length(schedules_selected_in_Swan_schedules) ) { 
        cat("Error: some schedules selected in Swan_schedules_combine.csv have no functions selected in Swan_function_list.csv", "\n") 
      }  
      
    } # end if Swan_function_list exists
    
    #   test if any species*countries are selected in Swan_schedules
    species_countries_selected <- Swan_schedules[, "B_sel"] == 1
    if (sum(species_countries_selected) == 0) { cat("Error! No species*countries selected in Swan_schedules.csv", "\n")  }
    
    #   ensure proper sorting of Swan_schedules.csv = sorting from highest Level to Lowest level  
    nlevel <- length(LevelCols)  # take into account that no of Levels might differ from the default 4 levels 
    
    if (nlevel == 3)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if (nlevel == 4)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"],Swan_schedules[,"Level4"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if (nlevel == 5)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"],Swan_schedules[,"Level5"], Swan_schedules[,"Level4"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if (nlevel == 6)  { 
      Swan_schedules                <- Swan_schedules[order(Swan_schedules[, "Schedule_code"], Swan_schedules[, "Species_nr"],Swan_schedules[,"Level6"],Swan_schedules[,"Level5"], Swan_schedules[,"Level4"], Swan_schedules[,"Level3"], Swan_schedules[,"Level2"], Swan_schedules[,"Level1"], decreasing=FALSE),]
    }
    if ((nlevel < 3) | (nlevel > 6)) { cat("Warning: ensure that the sorting of Swan_schedules_combine.csv is adequate", "\n") }
    
    Swan_schedules$CombiCol1	<- NULL
    Swan_schedules$CombiCol2	<- NULL
    Swan_schedules$CombiCol3	<- NULL
    write.table(Swan_schedules, file = "Swan_schedules_combine.csv",  dec=".",sep=";",row.names=FALSE) # properly sorted #13/02/2023 JAVI, DECIMAL SYMBOL
    
  } # end if Swan_schedules exists
  
  #Check for missing species-country combinations in species_countries.csv that are present in the schedule.csv #Javi 16/02/2023 new controls to prevent errors in rswan execution
  
  swan_schedule=read.table("Swan_schedules_combine.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
  Swan_schedules_selected <- Swan_schedules[Swan_schedules$B_sel == 1,]  
  speciescountry_code<-paste(species_countries_file$Species_nr,species_countries_file$Country_code)
  schedule_speciescountry <- paste(Swan_schedules_selected$Species_nr,Swan_schedules_selected$Level1)
  
  missing_combinations<-schedule_speciescountry[!(schedule_speciescountry %in% speciescountry_code) ]
  if(length(missing_combinations)!=0){
    missing_combinations <- gsub(" ","-",missing_combinations)
    cat("Error: The following species-country combination are selected in the schedule but not present in the species country file","\n" )
    cat("\n")
    cat(missing_combinations,"\n")
  } # END Check for missing species-country combinations # Javi 16/02/2023
  
  
} # E N D of function Check_steering_files
##########################################################################################################################################################
#END 20/04/2023 Javi




Run_Swan <- function(general_folder, working_folder,combine=0)  
  { 
  #combine=1
  ###################################################################################################
  # Aim: Main program Swan                                                                          #
  # Input: Steering files                                                                           #
  # Output: Output of functions which are called in function list per schedule                      #
  ###################################################################################################
  
  while (sink.number()  > 0 ) { sink() } 
  
  library(rtrim)
  
  setwd(general_folder) 
  
  if(combine ==0){#16/03/2023 Javi specifying the effect of combine argument
    comb<<-0
  ###########################################################################
  # read steering files                                                     #
  ###########################################################################
  
  # all schedules; selection of schedules occurs here
  Swan_schedules_name     <- "Swan_schedules.csv" 
  Swan_schedules          <<- read.table(Swan_schedules_name, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol    
  Swan_schedules_selected <- Swan_schedules[Swan_schedules$B_sel == 1,]                # delivers schedule file with selected schedules only
  schedules_selected_in_Swan_schedules <- as.character(unique(Swan_schedules[Swan_schedules$B_sel == 1,"Schedule_code"]))  
  
  LevelCols <- as.vector(grep("Level", names(Swan_schedules_selected)))
  # look up number of levels in Swan_schedules.csv (similar for all schedules)
  nlevels <- length(LevelCols)  # cannot be passed to called functions
  
  # function list, with selection of functions per schedule 
  Swan_function_list      <<- read.table("Swan_function_list.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
  
  # additional steering files 
  schedules_file          <<- read.table("Schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
  pools_file              <<- read.table("Pools.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol    
  countries_file          <<- read.table("Countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol     
  species_file            <<- read.table("Species.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol    
  species_countries       <<- read.table("Species_countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
  #species_countries$`Population size (geomean)`<<- gsub(",",".",species_countries$`Population size (geomean)`)#Javi 4/10/2021 there is an issue with the comma read as character vector
  }
  if(combine==1){#16/03/2023 Javi specifying the effect of combine argument
    comb<<-1
    # all schedules; selection of schedules occurs here
    Swan_schedules_name     <- "Swan_schedules_combine.csv" 
    Swan_schedules          <<- read.table(Swan_schedules_name, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol    
    Swan_schedules_selected <- Swan_schedules[Swan_schedules$B_sel == 1,]                # delivers schedule file with selected schedules only
    schedules_selected_in_Swan_schedules <- as.character(unique(Swan_schedules[Swan_schedules$B_sel == 1,"Schedule_code"]))  
    
    LevelCols <- as.vector(grep("Level", names(Swan_schedules_selected)))
    # look up number of levels in Swan_schedules.csv (similar for all schedules)
    nlevels <- length(LevelCols)  # cannot be passed to called functions
    
    # function list, with selection of functions per schedule 
    Swan_function_list      <<- read.table("Swan_function_list.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    
    # additional steering files 
    schedules_file          <<- read.table("Schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    pools_file              <<- read.table("Pools.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol    
    countries_file          <<- read.table("Countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol     
    species_file            <<- read.table("Species.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol    
    species_countries       <<- read.table("Species_countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    
    
    }#16/03/2023 END specifying the effect of combine argument
  ###########################################################################
  # switch to working folder, e.g PECBMS_2018 & open logfile                # 
  ###########################################################################
  working_folder <<- working_folder
  setwd(working_folder) 
  
  # start logfile
  sink("Logfile Swan.txt")
  cat("##################################","\n")
  cat(date(),"\n")  
  cat("\n")
  cat("Start Swan in R", "\n")  
  cat("\n")
  
  ###########################################################################
  # run schedules                                                           #
  ###########################################################################
  
  for (sch in 1:length(schedules_selected_in_Swan_schedules)) { # each schedule produces its own Swan-output (though in the same outputfolder)
    
    schedule_code        <<- schedules_selected_in_Swan_schedules[sch]  
    schedule_name        <<- as.character(schedules_file[schedules_file$Schedule_code == schedule_code, "Schedule_name"]) 
    Swan_schedule_active <<- Swan_schedules_selected[Swan_schedules_selected$Schedule_code == schedule_code,]  # << implies: pass value to called functions
    
    highest_year_of_schedule <<- 1950 # startvalue; will be adapted in function check_missing_years
    lowest_year_of_schedule  <<- 2050 # startvalue; will be adapted in function check_missing_years
    
    cat("Schedule applied:", schedule_code, schedule_name, "\n")
    cat("Working folder used:", working_folder, "\n")
    cat("\n")
    
    Swan_function_list_selected <- Swan_function_list[Swan_function_list$Schedule_code == schedule_code & Swan_function_list$B_sel == 1,]
    
    cat("The following functions in Swan_function_list will be executed:", "\n")
    for (f in 1:nrow(Swan_function_list_selected)) { 
      function_name         <- as.character(Swan_function_list_selected[f, "Function"])
      cat("function", function_name, "\n") 
    }
    
    ###########################################################################
    #   call selected functions from function list in the schedule            #
    ###########################################################################
    
    for (f in 1:nrow(Swan_function_list_selected)) { 
      function_name      <- as.character(Swan_function_list_selected[f, "Function"])
      numeric_arguments  <- Swan_function_list_selected[f, "Numeric_arguments"]
      arguments <- NULL
      arguments_cols       <- as.vector(grep("Argument", names(Swan_function_list_selected))) 
      number_arguments_cols <- length(arguments_cols)
      for (n in arguments_cols[1]:arguments_cols[number_arguments_cols]) { 
        if (!is.na(Swan_function_list_selected[f, n])) {
          if (numeric_arguments == TRUE)  {  arguments <- c(arguments, as.numeric(as.character(Swan_function_list_selected[f, n])))     }
          if (numeric_arguments == FALSE) {  arguments <- c(arguments, as.character(Swan_function_list_selected[f, n]))   }
        }
      }
      do.call(function_name, list(arguments))
      
    } # loop of functions within schedule
    
    cat("\n")
    cat("#############################################################################################################", "\n")
    cat("#############################################################################################################", "\n")
    cat("\n")
    
  } # end loop of schedules
  
  sink()
  
} # E N D of function Run_Swan
##########################################################################################################################################################








Select_all_available_datafiles <- function(no_arguments) { 
  ########################################################################################
  # Aim: Select all datafiles available in the working folder                            #
  # For simplicity's sake this is performed for ALL schedules at once                    #
  # Input: Datafiles in working_folder                                                   #  
  # Output: B_sel in Swan_schedules set at 1 if datafile appears to be present           #
  # Note: Swan_schedule already needs to have a record for the species & the stratum     #
  ########################################################################################
  
  cat("\n") 
  cat("@@@ function Select_all_available_datafiles", "\n")
  cat("\n") 
  
  Swan_schedules[, "B_sel"] <- 0 # all country files are deselected at the start
  no_files <- 1                   # vector with 1 indicating number of times of running loop only saved in R
  
  for (recordnr in 1: nrow(Swan_schedules)) { 
    Species       <- Swan_schedules[recordnr, "Species_nr"]
    Country       <- Swan_schedules[recordnr, "Level1"]
    datafile_name <- paste(Species, "_1_", Country, "_arg_output.csv",  sep="")
    if (file.exists(datafile_name)) { 
      Swan_schedules[recordnr, "B_sel"] <- 1  
      no_files <- no_files + 1
    } 
  }
  
  cat("B_sel in Swan_schedules set to 1 for ", no_files, "countries", "\n")   
  
  setwd(general_folder) 
  if(comb==1){#16/03/2023 Javi specifying the effect of combine argument
    write.table(Swan_schedules, file = "Swan_schedules_combine.csv",  dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
    Swan_schedules          <<- read.table("Swan_schedules_combine.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    setwd(working_folder) 
  }
  
  if (comb==0){
  write.table(Swan_schedules, file = "Swan_schedules.csv",  dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  
  Swan_schedules          <<- read.table("Swan_schedules.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
  setwd(working_folder) 
  }#16/03/2023 END specifying the effect of combine argument
  cat("\n")
  
} # E N D of function Select_all_available_datafiles
##########################################################################################################################################################



Check_existence_datafiles <- function(selection_of_files) { # if TRUE then B_sel set at 1 if data are available
  ########################################################################################
  # Aim: Test if all RTRIM outputfiles of countries selected in Swan_schedules.csv exist #
  # Input: List of RTRIM outputfiles of countries (= level1)                             #  
  # Output: Errors in logfile                                                            #
  # Note: Schedule is not in the name of the RTRIM outputfiles                           # 
  ########################################################################################
  
  cat("\n") 
  cat("@@@ function Check_existence_datafiles executed", "\n")
  cat("\n") 
  
  if (nrow(Swan_schedule_active) == 0) { cat("Error! No selection of species and strata in", schedule_code,"\n") }  
  
  if (nrow(Swan_schedule_active) > 0) { 
    
    ###########################################################################
    #   test if datafiles are missing in working_folder                       #
    ###########################################################################
    
    cat("test how many species_countries selected in Swan_schedules are not available in the working_folder", "\n")
    
    no_of_datafiles_found_in_folder <- 0
    listSpeciesStratumCombinations_in_schedule <- NULL
    
    for (r in 1:nrow(Swan_schedule_active)) { # Note: loop follows the records in Swan_schedule which have been selected (B_sel = 1)
      
      Species     <- Swan_schedule_active[r, "Species_nr"]
      Country     <- Swan_schedule_active[r, "Level1"]
      
      RTRIM_timetotals   <-  paste(Species, "_1_", Country, "_indices_TT.csv",  sep="")  
      if (file.exists(RTRIM_timetotals)  == FALSE)   { cat("Warning: datafile", RTRIM_timetotals, "missing in folder", "\n")   }
      if (file.exists(RTRIM_timetotals)  == TRUE)    { no_of_datafiles_found_in_folder <- no_of_datafiles_found_in_folder + 1  }
      RTRIM_ocv          <-  paste(Species, "_1_", Country, "_ocv.csv",  sep="")   
      if (file.exists(RTRIM_ocv)  == FALSE)          { cat("Warning: datafile", RTRIM_ocv, "missing in folder", "\n")          }
      if (file.exists(RTRIM_ocv)  == TRUE)           { no_of_datafiles_found_in_folder <- no_of_datafiles_found_in_folder + 1  }
      RTRIM_arguments    <-  paste(Species, "_1_", Country, "_arg_output.csv",  sep="")
      if (file.exists(RTRIM_arguments)  == FALSE)    { cat("Warning: datafile", RTRIM_arguments, "missing in folder", "\n")    }
      if (file.exists(RTRIM_arguments)  == TRUE)     { 
        no_of_datafiles_found_in_folder <- no_of_datafiles_found_in_folder + 1         
        #       put name on the list; for simplicy only arg_output files incorporated in the list 
        listSpeciesStratumCombinations_in_schedule <- c(listSpeciesStratumCombinations_in_schedule, RTRIM_arguments)
      }
      
    } # test inclusion of species_countries in Swan_schedules
    
    cat("number of datafiles selected in working folder:", 3 * nrow(Swan_schedule_active), "\n") 
    cat("number of datafiles found in working folder   :", no_of_datafiles_found_in_folder, "\n")
    cat("number of datafiles selected in schedule, but not found in working_folder:", (3 * nrow(Swan_schedule_active)) - no_of_datafiles_found_in_folder, "\n")
    cat("\n") 
    
    ###########################################################################
    #   test if datafiles are perhaps missing in Swan_schedules               #
    ###########################################################################
    
    cat("test how many species_countries available in the working_folder are not included in Swan_schedules (even if B_sel == 0)", "\n")
    #   only arg_output files are incorporated in this test; they serve here as proxy for all species_country files 
    
    no_of_datafiles_not_in_schedule <- 0
    listSpeciesStratumCombinations_in_folder <- NULL
    
    listSpeciesStratumCombinations_in_folder <- dir(working_folder, pattern = "arg_output.csv")
    #   filter only country files
    #   name of country datafile = <euring>-<1>-<stratum>.csv                thus contains two hyphens
    #   name of pool datafile    = <schedule>-<euring>-<1>-<stratum> .csv    thus contains three hyphens
    
    for (s_st in 1: length(listSpeciesStratumCombinations_in_folder)) { 
      #     remove tekst from filename until next hyphen
      remaining_text <- gsub("^[a-zA-Z0-9]*_", "", listSpeciesStratumCombinations_in_folder[s_st])   
      remaining_text <- gsub("^[a-zA-Z0-9]*_", "", remaining_text)
      remaining_text <- gsub("^[a-zA-Z0-9]*_", "", remaining_text)
      if (remaining_text == "arg_output.csv") { # datafile found belonging to a country (= having two hyphens only)
        presence_in_schedule <- is.element(listSpeciesStratumCombinations_in_folder[s_st], listSpeciesStratumCombinations_in_schedule)
        if (presence_in_schedule == FALSE) { 
          cat("Warning: datafile", listSpeciesStratumCombinations_in_folder[s_st], "missing in Swan_schedule.csv", "\n")          
          no_of_datafiles_not_in_schedule <- no_of_datafiles_not_in_schedule + 1
        }
      }
    }
    
    cat("number of datafiles present in working folder, but not in Swan_schedule:", no_of_datafiles_not_in_schedule, "\n")
    
  } # end if any selected records in schedule
  
} # E N D of function Check_existence_datafiles
##########################################################################################################################################################








Check_content_datafiles <- function(no_arguments) { 
  ###########################################################################
  # Aim: Test content of all RTRIM outputfiles of countries                 # 
  # Checks are similar as in TOIP                                           #
  # Input: RTRIM output files of countries                                  #  
  # Output: Errors in logfile                                               #
  # Note: Schedule is not in the name of the RTRIM outputfiles              # 
  ###########################################################################
  
  cat("\n") 
  cat("@@@ function Check_content_datafiles executed", "\n")
  cat("\n") 
  
  min_year <- 2100 # starting value; will be adapted in this function 
  max_year <- 0    # starting value; will be adapted in this function 
  
  for (r in 1:nrow(Swan_schedule_active)) { # Note: loop follows the records in Swan_schedule
    
    Species     <- Swan_schedule_active[r, "Species_nr"]
    Country     <- Swan_schedule_active[r, "Level1"]
    
    #   start to check content indices_TT outputfile
    RTRIM_timetotals   <- paste(Species, "_1_", Country, "_indices_TT.csv",  sep="")  
    #print(paste0("Estoy intentando abrir ",RTRIM_timetotals))
    TimetotalsStratum  <- read.table(RTRIM_timetotals, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    
    #   expecting 9 columns: Year TT_model TT_model_SE TT_imputed TT_imputed_SE Index_model Index_model_SE Index_imputed Index_imputed_SE
    if (ncol(TimetotalsStratum) != 9)              { cat("Warning: datafile", RTRIM_timetotals, "has an unexpected number of columns", "\n")      }
    
    ###########################################################################
    #   check for suspicious indices and time totals                          #
    ###########################################################################
    
    low_indices_present <- TimetotalsStratum[, "Index_imputed"] < 0.005 # note: index 0.005 stands for index 0.5 if base year index is 100
    if (sum(low_indices_present > 0))              { cat("Warning: datafile", RTRIM_timetotals, "contains indices < 0.5", "\n")                   }         
    
    high_indices_present <- TimetotalsStratum[, "Index_imputed"] > 10 # note: index 10 stands for index 1000 if base year index is 100
    if (sum(high_indices_present > 0))             { cat("Warning: datafile", RTRIM_timetotals, "contains indices > 1000", "\n")                  }         
    
    low_timetotals_present <- TimetotalsStratum[, "TT_imputed"] < 1 
    if (sum(low_timetotals_present > 0))           { cat("Warning: datafile", RTRIM_timetotals, "contains time totals < 0.5", "\n")               }         
    
    high_timetotals_present <- TimetotalsStratum[, "TT_imputed"] > 10000 
    if (sum(high_timetotals_present > 0))          { cat("Warning: datafile", RTRIM_timetotals, "contains time totals > 1000", "\n")              }         
    
    low_years_present <- TimetotalsStratum[, "Year"] < 1990
    if (sum(low_years_present > 0))                { cat("Warning: datafile", RTRIM_timetotals, "contains years before 1990", "\n")               }         
    
    current_year <- as.numeric(format(Sys.time(), "%Y"))
    high_years_present <- TimetotalsStratum[, "Year"] > current_year 
    if (sum(high_years_present > 0))               { cat("Warning: datafile", RTRIM_timetotals, "contains years after", current_year, "\n")       }         
    
    ###########################################################################
    #   basetime checks                                                       #
    ###########################################################################
    
    index_base_year     <- TimetotalsStratum[TimetotalsStratum$Index_imputed == 1, "Year"]
    
    
    if (length(index_base_year) == 0 ) { cat("Warning: datafile", RTRIM_timetotals, "contains no basetime indices", "\n")       }  
    
    if (length(index_base_year)>= 2 ) {
      cat("Warning: datafile", RTRIM_timetotals, "contains more years with indices = 1 than the base year", "\n")# 13/02/2023 If for some reason there are more than a year with an index value of 1 it will print a warning and will select the first year among them
      index_base_year<-index_base_year[1] 
    }
    
    if (length(index_base_year) >  0 ) { 
      index_first_year    <- TimetotalsStratum[1, "Year"]
      #print(paste0 ("el index base year es: ", index_base_year," y el index first year es: ",index_first_year))
      if (index_base_year != index_first_year)       { cat("Warning: datafile", RTRIM_timetotals, "indices basetime is not first year", "\n")       }         
      
      base_years_present  <- TimetotalsStratum[(TimetotalsStratum$Index_imputed == 1) & (TimetotalsStratum$Index_imputed_SE == 0), "Year"]
      if (length(base_years_present) > 1)            { cat("Warning: Datafile", RTRIM_timetotals, "contains >1 year with index=100 and se=0", "\n") }     
    }
    
    ###########################################################################
    #   assess span of years of all countries together                        # 
    ###########################################################################
    
    years_present       <- TimetotalsStratum[, "Year"]
    
    max_year_stratum <- max(years_present) 
    if (max_year_stratum > max_year) { max_year <<-  max_year_stratum }    # max_year is for further use in another function
    
    min_year_stratum <- min(years_present) 
    if (min_year_stratum < min_year) { min_year <<-  min_year_stratum }    # max_year is for further use in another function
    
    ###########################################################################
    #   check content of ocv outputfile                                       # Arco 18.12.2020
    ###########################################################################
    
    RTRIM_ocv          <-  paste(Species, "_1_", Country, "_ocv.csv",  sep="") 
    #print(paste0("Estoy intentando abrir ",RTRIM_ocv))
    TimeocvStratum  <- read.table(RTRIM_ocv, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
    
    #   test if ocv file is symmetrical
    nrows    <- nrow(TimeocvStratum)
    ncolumns <- ncol(TimeocvStratum)
    if (nrows != ncolumns) { cat("Error: ocv file is asymmetrical", "\n")  }
    
    #   test if ocv and tt are about the same species and stratum = they have (almost) similar year variances
    # Line: "year year_variance_in_tt  <- as.numeric(TimetotalsStratum[c_ocv, "TT_imputed_SE"]) * as.numeric(TimetotalsStratum
    #  [c_ocv, "TT_imputed_SE"])  # Arco 17.12.2020" was added by Arco
    # = code solving integer overflow, because of the extreme large numbers involved. 
    for (r_ocv in 1:nrows) {
      for (c_ocv in 1:ncolumns) {
        if (r_ocv == c_ocv) {
          year_variance_in_ocv <- TimeocvStratum [r_ocv, c_ocv] 
          year_variance_in_tt  <- as.numeric(TimetotalsStratum[c_ocv, "TT_imputed_SE"]) * as.numeric(TimetotalsStratum[c_ocv, "TT_imputed_SE"])  # Arco 17-12-2020
          mutual_deviation_in_se <- ((year_variance_in_ocv - year_variance_in_tt) / year_variance_in_ocv) * 100
          if (mutual_deviation_in_se > 25) {
            cat("Warning:", RTRIM_timetotals, "and", RTRIM_ocv, "have > 25% difference in variances of time totals", "\n")
          }
        }
      }
    }
    
    
    ###########################################################################
    #   check content of arg_output file                                      #
    ###########################################################################
    
    RTRIM_arguments    <- paste(Species, "_1_", Country, "_arg_output.csv",  sep="")
    
    #   cat("RTRIM_arguments", RTRIM_arguments, "\n") # only for testing
    #print(paste0("Estoy intentando abrir ",RTRIM_arguments))
    TimeargStratum     <- read.table(RTRIM_arguments, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol   
    
    #   Expecting 18 columns: N_sites; N_time_values; N_observed_zero_counts; N_observed_positive_counts; N_missing_counts; N_counts_total; Base_year_first_year; Base_year_last_year; Changepoints; 
    #   Overdispersion; Serial_correlation; Slope_imputed_mul; Slope_imputed_mul_SE; Slope_imputed_classification; Year_from; Slope_from_imputed_mul; Slope_from_imputed_mul_SE; Slope_from_imputed_classification  
    # slightly adapted 29-1-2019 checked by ES 15-1-2021
    if (ncol(TimeargStratum) != 18)              { cat("Warning: datafile", RTRIM_arguments, "has an unexpected number of columns", "\n")      }
    
    # If the problem is the "NA" in column Overdispersion in _agr_output.csv (Correction by Javier 31.12021)
    #if(is.na(TimeargStratum[,which(names(TimeargStratum)=="Overdispersion")])) {TimeargStratum$Overdispersion<-1}
    #Javi04/11/2021 removed cause data inputs do not have this column in the analyses of 2019
    
    ##### ES: Dutch files are done in different way and a column "Overdispersion was missing -> if this happens, 
    # we YOU SHOULD ASK FOR THE EXISTENCE OF OVERDISPERSION (USING COMMAND "WHICH"). 
    #IF IT IS ABSENT CREATE COLUMN: RTRIM_arguments$Overdispersion <- 1 (or a similar code)
    # forward Arco's email from 11/12
    # I want to have the line here with slash for case of need x this is second control for missing columns in _output.csv (first is in ONLT)
    # => The line will be behind hashtag
    ###########################################################################
    #   test for multiple runs in one (old) TRIM output                       #
    #   TOIP then will generate more rows in the argfile                      #     
    ###########################################################################
    
    if (nrow(TimeargStratum) != 1)               { cat("Warning: datafile", RTRIM_arguments, "has more than one row, probably because multiple (old) TRIM runs were available in the outfile", "\n")      }
    
    ###########################################################################
    #   check for suspicious time totals (counts)                             #
    ###########################################################################
    
    n_zero_counts_present        <- TimeargStratum[, "N_observed_zero_counts"] < 1 
    if (n_zero_counts_present == TRUE) { cat("Warning: the RTRIM analysis to produce", RTRIM_arguments, "contained no zero counts (perhaps RSWAN has been applied?)", "\n")       }         
    
    n_missing_counts_present     <- TimeargStratum[, "N_missing_counts"] < 1 
    if (n_missing_counts_present == TRUE) { cat("Warning: the RTRIM analysis to produce", RTRIM_arguments, "contained no missing counts (perhaps RSWAN has been applied?)", "\n") }         
    
    #overdispersion_applied       <- TimeargStratum[, "Overdispersion"] != 1
    #if (is.na(overdispersion_applied)) { 
    #  cat("Warning: the RTRIM analysis to produce", RTRIM_arguments, "applied no overdispersion", "\n")          
    #} else {
    #  if (overdispersion_applied == FALSE) { cat("Warning: the RTRIM analysis to produce", RTRIM_arguments, "applied no overdispersion", "\n")   }       
    #}
    
    #serial_correlation_applied   <- TimeargStratum[, "Serial_correlation"] != 0
    #if (is.na(serial_correlation_applied == FALSE)) { 
    #  cat("Warning: the RTRIM analysis to produce", RTRIM_arguments, "applied no serial correlation", "\n")     
    #} else { 
    #  if (serial_correlation_applied == FALSE) { cat("Warning: the RTRIM analysis to produce", RTRIM_arguments, "applied no serial correlation", "\n") }         
    #}
    
    ###########################################################################
    #   check if linear trend model with no changepoints has been used        #
    ###########################################################################
    
    changepoints_applied         <- as.character(TimeargStratum[, "Changepoints"])
    number_of_changepoints       <- unlist(lapply(strsplit(changepoints_applied, split = "-"), length))  # count the number of elements in the string
    if (number_of_changepoints + 1 == nrows) { cat("All changepoints were used in the model to produce", RTRIM_arguments, "\n") }         
    if (number_of_changepoints <= 2)  { cat("Warning: the linear trend model to produce", RTRIM_arguments, "used too few changepoints", "\n") }         
    
    ###########################################################################
    #   check for suspicious slopes                                           #
    ###########################################################################
    
    low_slope_present             <- TimeargStratum[, "Slope_imputed_mul"] < 0.6
    if (low_slope_present)        { cat("Warning: slope estimate in ", RTRIM_arguments, "< 0.6", "\n")              }         
    
    high_slope_present            <- TimeargStratum[, "Slope_imputed_mul"] > 1.5
    if (high_slope_present)       { cat("Warning: slope estimate in ", RTRIM_arguments, "> 1.5", "\n")              }         
    
    low_slope_from_present        <- TimeargStratum[, "Slope_from_imputed_mul"] < 0.6
    if (!is.na(low_slope_from_present)) {  
      if (low_slope_from_present)   { 
        cat("Warning: slope_from estimate in ", RTRIM_arguments, "< 0.6", "\n")   
      }
    }         
    
    high_slope_from_present       <- TimeargStratum[, "Slope_from_imputed_mul"] > 1.5
    if (!is.na(high_slope_from_present)) {  
      if (high_slope_from_present)  { 
        cat("Warning: slope_from estimate in ", RTRIM_arguments, "> 1.5", "\n")     
      }
    }         
    
    low_slope_error_present       <- TimeargStratum[, "Slope_imputed_mul_SE"] < 0.001
    if (low_slope_error_present)  { cat("Warning: slope error estimate in ", RTRIM_arguments, "< 0.001", "\n")      }         
    
    high_slope_error_present      <- TimeargStratum[, "Slope_imputed_mul_SE"] > 0.5
    if (high_slope_error_present) { cat("Warning: slope error estimate in ", RTRIM_arguments, "> 0.5", "\n")        }         
    
  } # end nrow(Swan_schedules_sel)
  
} # E N D of function Check_content_datafiles
##########################################################################################################################################################





Check_missing_years_in_level2_pools <- function(no_arguments) {  
  ####################################################################################
  # Aim: Test if missing years are expected in the level2 pools                      # 
  # Note: Errors are not checked in higher (>2) level pools                          #
  # Input: indices_TT.csv files of the countries that are part of the pool           # 
  # Output: Errors in Logfile                                                        # 
  ####################################################################################
  
  cat("\n")
  cat("@@@ function Check_missing_years_in_level2_pools executed", "\n") 
  cat("\n") 
  
  LevelCols <- as.vector(grep("Level", names(Swan_schedules)))
  # look up number of levels in Swan_schedules.csv (similar for all schedules)
  
  Pool_column                     <- paste("Level", 2, sep ="")
  Pool_column_before              <- paste("Level", 1, sep ="")   # country column
  
  all_species <- unique(Swan_schedule_active[Swan_schedule_active$B_sel == 1, "Species_nr"])
  nspecies    <- length(all_species)
  
  for (s in 1:nspecies) { 
    
    pools     <- unique(Swan_schedule_active[Swan_schedule_active$Species_nr == all_species[s], Pool_column])
    npools    <- length(pools)
    
    for (p in 1:npools) { 
      
      #     strata are pools of lower level
      strata       <- unique(Swan_schedule_active[Swan_schedule_active$Species_nr == all_species[s] & Swan_schedule_active[, Pool_column] == pools[p], Pool_column_before])
      nstrata      <- length(strata)
      
      ###########################################################################
      #     determine number of years within pool                             #
      ###########################################################################
      
      allyear <- NULL 
      
      for (st in 1:nstrata) { 
        
        filenameTimeTotals_out <- paste(all_species[s], "_",  1, "_",  strata[st], "_indices_TT.csv", sep="")  # level1 = country files input; no schedule_code in the name 
        TimetotalsStratum  <- read.table(filenameTimeTotals_out, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol  
        stratum_year       <- unique(TimetotalsStratum[, "Year"])
        allyear            <- c(allyear, stratum_year) 
        
      } # end loop strata
      
      uyear    <- sort(unique(allyear))	
      nyear    <- length(uyear)
      
      maxyear  <- as.numeric(max(uyear)) # maxyear of pool
      minyear  <- as.numeric(min(uyear)) # minyear of pool
      
      if ((maxyear-minyear+1) != nyear) { cat("Error: pool", pools[p], "of species", all_species[s], "has missing years in all lower level strata taken together", "\n") }
      if ((maxyear-minyear+1) == nyear) { cat("Ok: pool", pools[p], "of species", all_species[s], "has no missing years in all lower level strata taken together", "\n") }
      
    } # end loop pools
    
  } # end loop species
  
} # E N D of function Check_missing_years_in_level2_pools
##########################################################################################################################################################









Import_weights_into_Swan_schedules <- function(no_arguments) { 
  ############################################################################
  # Aim: Calculate weights and import these into Swan_schedules.csv          #
  # Input: Time totals per species & estimates of population size / country  #
  # Output: updated weights in Swan_schedules.csv                            #
  ############################################################################
  
  cat("\n") 
  cat("@@@ function Import_weights_into_Swan_schedules executed", "\n")
  cat("\n") 
  
  for (r in 1:nrow(Swan_schedule_active)) { 
    
    Species     <- Swan_schedule_active[r, "Species_nr"]
    
    Country     <- Swan_schedule_active[r, "Level1"]
    
    
    RTRIM_timetotals   <- paste(Species, "_1_", Country, "_indices_TT.csv",  sep="") 
    TimetotalsStratum  <- read.table(RTRIM_timetotals, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
    
    ###########################################################################
    #   check consisteny time totals with population sizes                    #
    ###########################################################################
    #   population sizes in species_country table (with estimates from e.g. 1999-2001) 
    
    species_country_mean_popsize     <- as.numeric(species_countries[species_countries$Species_nr == Species & species_countries$Country_code == Country, "Population size (geomean)"])
    species_country_max_popsize      <- as.numeric(species_countries[species_countries$Species_nr == Species & species_countries$Country_code == Country, "Population size_max"])
    species_country_popsize_first_yr <- species_countries[species_countries$Species_nr == Species & species_countries$Country_code == Country, "Weight_Year_First"]
    species_country_popsize_last_yr  <- species_countries[species_countries$Species_nr == Species & species_countries$Country_code == Country, "Weight_Year_Last"]
    if (Country== 64 | Country==65 |Country == 69| Country ==70 | Country==71 |Country==72){ #16-3-2023 Javi, added this lane in Feb when building the combine function to dive population of germany E,W and sweden to the half
      species_country_mean_popsize<- species_country_mean_popsize/2
      species_country_max_popsize<- species_country_max_popsize/2
    }
    #   test for mismatch of years 
    species_country_timetotal <- mean(TimetotalsStratum[TimetotalsStratum$Year >= species_country_popsize_first_yr & TimetotalsStratum$Year <= species_country_popsize_last_yr,"TT_imputed"] )
    if (is.na(species_country_timetotal) == TRUE) { 
      cat("Error! Mismatch of years in data and popsize estimate for species", Species, "in country", Country, "\n")
    } 
    
    if (is.na(species_country_timetotal) == FALSE) { 
      ProcPopSizeTimeTotals          <- (species_country_timetotal/species_country_mean_popsize)*100
      if (ProcPopSizeTimeTotals > 95)                              { cat("Warning: datafile", RTRIM_timetotals, "time totals > 95% of mean population size", "\n")    }         
      if (species_country_timetotal > species_country_max_popsize) { #Javi 04/10/2021 I've modifyed the function to place 1 as weight when the TT are hieguer than the max population 
        cat("Warning: datafile", RTRIM_timetotals, "time totals > max population size", "\n")  
         Weight=1}
      else{
        Weight          <- species_country_mean_popsize/species_country_timetotal # inverse of ProcPopSizeTimeTotals
        }      
      
      ###########################################################################   
      #     calculate weight & transport to Swan_schedules.csv col "Weight"     #
      ###########################################################################
      
      
      Swan_schedules[Swan_schedules$Schedule_code == schedule_code & Swan_schedules$Species == Species & Swan_schedules$Level1 == Country, "Weight"] <- Weight
      
      #     weights need also to be adapted in the active schedule, else new weights will not be used in Rtrimcovin:
      Swan_schedule_active[Swan_schedule_active$Schedule_code == schedule_code & Swan_schedule_active$Species == Species & Swan_schedule_active$Level1 == Country, "Weight"] <- Weight
      #print(paste0("the weight stored is ",Swan_schedule_active[Swan_schedule_active$Schedule_code == schedule_code & Swan_schedule_active$Species == Species & Swan_schedule_active$Level1 == Country, "Weight"]))
      }
    
  } # end nrow(Swan_schedules_sel)
  
  setwd(general_folder)
  if(comb==1){#16/03/2023 Javi specifying the effect of combine argument
    write.table(Swan_schedules, file = "Swan_schedules_combine.csv",  dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  }
  if(comb==0){
    write.table(Swan_schedules, file = "Swan_schedules.csv",  dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  }#16/03/2023 END specifying the effect of combine argument
  cat("Weights calculated and written in column Weight in Swan_schedules.csv", "\n")
  setwd(working_folder) 
  Swan_schedule_active <<-  Swan_schedule_active  # ARCO 3-10-2021
  
} # E N D of function Import_weights_into_Swan_schedules
##########################################################################################################################################################




Truncate_datafiles <- function(no_arguments) { 
  ###########################################################################
  # Aim: Truncate RTRIM outputfiles which are active in SWAN_schedules.csv  #
  # Input: RTRIM shell files tt & ocv                                       #
  # Output: truncated RTRIM shell files tt & ocv; argfiles remain unchanged #
  ###########################################################################
  
  cat("\n") 
  cat("@@@ function Truncate_datafiles", "\n")
  cat("\n") 
  
  Swan_active_records <- Swan_schedule_active[Swan_schedule_active$B_sel == 1, ]
  nrecords            <- nrow(Swan_active_records)
  truncated_files     <- 0
  
  for (st in 1:nrecords) { 
    
    species_nr <- Swan_active_records[st, "Species_nr"]
    stratum_nr <- Swan_active_records[st, "Level1"]
    
    #   look up first and last year in Species_Countries file
    year_first <- species_countries[species_countries$Species_nr == species_nr & species_countries$Country_code == stratum_nr, "Year_first"] 
    year_last  <- species_countries[species_countries$Species_nr == species_nr & species_countries$Country_code == stratum_nr, "Year_last"] 
    
    #   RTRIM time totals
    filenameTimeTotals_out <- paste(species_nr, "_",  1, "_",  stratum_nr, "_indices_TT.csv", sep="")  # level1 = country files input; no schedule_code in the name 
    TimetotalsStratum  <- read.table(filenameTimeTotals_out, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
    nrows_before <- nrow(TimetotalsStratum)
    range_years_trimoutput <- TimetotalsStratum[1, "Year"] : TimetotalsStratum[nrow(TimetotalsStratum), "Year"] # information needed for truncation of ocvfiles
    #   remove rows:
    TimetotalsStratum <- TimetotalsStratum[TimetotalsStratum$Year >= year_first & TimetotalsStratum$Year <= year_last,]
    nrows_after <- nrow(TimetotalsStratum)
    
    if (nrows_before > nrows_after) { # truncation has occurred
      cat("Datafile", filenameTimeTotals_out, "has been truncated by removing", nrows_before - nrows_after, "years", "\n") 
      write.table(TimetotalsStratum, file=filenameTimeTotals_out, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
      
      #     also truncate RTRIM ocv file:
      filenameCovmatrix <- paste(species_nr, "_",  1, "_",  stratum_nr, "_ocv.csv", sep="")  # level1 = country file input
      ocvStratum <- read.table(filenameCovmatrix, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
      #     ocv file does not have a column with calendar year, so use year information from range_years_trimoutput instead
      #     remove both rows and columns:
      row_colums_to_keep <- (range_years_trimoutput >= year_first) & (range_years_trimoutput <= year_last)
      ocvStratum <- ocvStratum[row_colums_to_keep,row_colums_to_keep]
      write.table(ocvStratum, file=filenameCovmatrix, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
      
      truncated_files <- truncated_files + 1
    }
    
  } # end loop strata
  
  cat("Number of truncated files:", truncated_files, "\n")  
  cat("\n")
  
} # E N D of function Truncate_datafiles
##########################################################################################################################################################






Strip_ocv_files <- function(ocvfiles_to_strip) { 
  ###########################################################################
  # Aim: Set year-year covariances at 0 to treat estimation difficulties    #
  # Input: Ocv files mentioned as arguments in Swan_function list           #
  # Output: Ocv files with all covariances having value 0                   #
  ###########################################################################
  
  cat("\n") 
  cat("@@@ function Strip_ocv_files", "\n")
  cat("\n") 
  
  for (i in 1:length(ocvfiles_to_strip)) {
    
    if (!is.na(ocvfiles_to_strip[i])) { 
      ocvtest <- read.table(ocvfiles_to_strip[i],header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol  # read ocv
      no_rows    <- nrow(ocvtest)
      no_cols    <- ncol(ocvtest)
      for (r_ocv in 1:no_rows) {
        for (c_ocv in 1:no_cols) {
          if (r_ocv != c_ocv) ocvtest[r_ocv, c_ocv] <- 0 # set all covariances at zero
        }
      }
    } 
    
    ###########################################################################   
    #   overwrite original ocv file                                           #
    ###########################################################################   
    write.table(ocvtest, file = ocvfiles_to_strip[i], dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
    cat("ocv file stripped:", ocvfiles_to_strip[i], "\n")
    
  } 
  
} # E N D of function Strip_ocv_files
##########################################################################################################################################################






Make_empty_arg_output_file <- function(no_arguments) {
  ################################################################################
  # Aim: this function creates an empty dataframe for slopes spanning the entire #
  # study period and a subperiod for a species-stratum combination.              # 
  # the arg_output_file is not needed for SWAN, but it makes the output similar  #
  # to output of RTRIM shell                                                     # 
  ################################################################################
  
  cat("\n") 
  cat("@@@ function Make_empty_arg_output_file", "\n")
  cat("\n") 
  
  arg_output_file <<- data.frame(N_sites = integer(1),        # Number of unique sites.
                                 N_time_values = integer(1),                      # Number of unique years.
                                 N_observed_zero_counts = integer(1),             # Number of zero counts.
                                 N_observed_positive_counts = integer(1),         # Number of positive counts.
                                 N_missing_counts = integer(1),                   # Number of missing counts.
                                 N_counts_total = integer(1),                     # Total number of counts.
                                 Base_year_first_year = numeric(1),
                                 Base_year_last_year = numeric(1),
                                 # Calendar year used as base year for indices. 
                                 # If Base_year_first_year equals Base_year_last_year a single year is used as base year.
                                 # If Base_year_first_year < Base_year_last_year, a period is used as base time.
                                 # In the latter case, Base_year_first_year is the first year of the period.
                                 
                                 Changepoints = numeric(1),                        # Changepoints used.
                                 Overdispersion = numeric(1),                      # Estimated overdispersion.
                                 Serial_correlation = numeric(1),                  # Estimated serial correlation.
                                 Slope_imputed_mul = numeric(1),                   # Multiplicative imputed slope for the entire period.
                                 Slope_imputed_mul_SE = numeric(1),                # Standard error of multiplicative imputed slope for the entire period.
                                 Slope_imputed_classification = character(1),      # Trend classification of multiplicative imputed slope for the entire period.
                                 Year_from = numeric(1),                           # First year of the subperiod from which a slope has been calculated. 
                                 # Last year of the subperiod is equal to the last year of the entire period.
                                 # Note that it is assumed that the first year of the subperiod is closer to the present than the first year of the entire period.
                                 
                                 Slope_from_imputed_mul = numeric(1),              # Multiplicative imputed slope for the subperiod.
                                 Slope_from_imputed_mul_SE = numeric(1),           # Standard error of multiplicative imputed slope for the subperiod.
                                 Slope_from_imputed_classification = character(1)  # Trend classification of multiplicative imputed slope for the subperiod.
                                 
  )
  
} # E N D of function Make_empty_arg_output_file
##########################################################################################################################################################




Rtrimcovin <- function(base_year) {  #Javi 21.09.2021
  ###############################################################################
  # Aim: Run RTRIM trimcovin for all levels                                     #
  # Input level 2  : outputfiles of countries (no schedule_code in the name)    # 
  # Input level > 2: time totals and icv files (schedule_code in the name)      # 
  # Output: RTRIM RData output (schedule in the name)                           #
  # Output: RTRIM TT and ocv files output (schedule in de name)                 #
  ###############################################################################
  
  cat("\n")
  cat("@@@ function Rtrimcovin executed", "\n") 
  cat("\n") 
  
  #######################################################################################################
  
  LevelCols <- as.vector(grep("Level", names(Swan_schedules)))
  # look up number of levels in Swan_schedules.csv (similar for all schedules)
  nlevels <- length(LevelCols)  
  
  for (level in 2:nlevels) { 
    
    Pool_column                     <- paste("Level", level, sep ="") 
    Pool_column_before              <- paste("Level", level - 1, sep ="") 
    
    all_species <- unique(Swan_schedule_active[Swan_schedule_active$B_sel == 1, "Species_nr"])
    nspecies    <- length(all_species)
    
    for (s in 1:nspecies) { 
      
      pools     <- unique(Swan_schedule_active[Swan_schedule_active$Species_nr == all_species[s], Pool_column])
      npools    <- length(pools)
      
      for (p in 1:npools) { 
        
        #       strata as pools of lower level
        strata       <- unique(Swan_schedule_active[Swan_schedule_active$Species_nr == all_species[s] & Swan_schedule_active[, Pool_column] == pools[p], Pool_column_before])
        nstrata      <- length(strata)
        
        ###########################################################################
        #       determine number of years within pool                             #
        ###########################################################################
        
        allyear <- NULL 
        for (st in 1:nstrata) { # jaren opzoeken
          
          if (level == 2) { 
            filenameTimeTotals_out <- paste(all_species[s], "_",  level-1, "_",  strata[st], "_indices_TT.csv", sep="")  # level1 = country file input 
            
          } 
          if (level > 2) { 
            filenameTimeTotals_out <- paste(schedule_code, "_", all_species[s], "_",  level-1, "_",  strata[st], "_indices_TT.csv", sep="")  # pool file input, with schedule_code in the namr
          }
          #print(paste0("Leyendo ",filenameTimeTotals_out))
          TimetotalsStratum  <- read.table(filenameTimeTotals_out, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
          stratum_year       <- unique(TimetotalsStratum[, "Year"])
          allyear            <- c(allyear, stratum_year) 
          
        } # end loop strata
        
        uyear    <<- sort(unique(allyear))
        nyear    <<- length(uyear)
        
        #       yearfile, containing all years of all strata
        yearfile <- array(NA,dim=c(nyear, 1)) 
        colnames(yearfile) <- c("Year")
        yearfile [, "Year"] <- uyear
        
        ###########################################################################
        #       combine time total files of all strata of pool together           #
        ###########################################################################
        
        TimetotalsAllStrata <- NULL
        
        for (st in 1:nstrata) {
          if (level == 2) { 
            filenameTimeTotals_out <- paste(all_species[s], "_",  level-1, "_",  strata[st], "_indices_TT.csv", sep="")  # level1 = country file input 
          } 
          if (level > 2) { 
            filenameTimeTotals_out <- paste(schedule_code, "_", all_species[s], "_",  level-1, "_",  strata[st], "_indices_TT.csv", sep="")  # pool file input, with schedule in the name 
          }
          #print(paste0("Leyendo ",filenameTimeTotals_out))
          TimetotalsStratum        <- read.table(filenameTimeTotals_out, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
          #         include possible missing years 
          TimetotalsStratum        <- merge(yearfile, TimetotalsStratum[, c("Year", "TT_imputed")], by = c("Year"), all = TRUE)
          
          TimetotalsStratum$site   <- strata[st] # stratumcode, could be country or pool
          
          #         fill in column count: 
          TimetotalsStratum$count  <- TimetotalsStratum$TT_imputed
          TimetotalsStratum$TT_imputed <- NULL
          
          TimetotalsStratum$weight <- 1  # default, plus required if level > 2 and no weighting in RTRIM needed
          
          if (level == 2) { # weight derived from Swan_schedules.csv
            weight_from_Swan_schedule <- as.numeric(Swan_schedule_active[Swan_schedule_active$Species == all_species[s] & Swan_schedule_active$Level1 == strata[st], "Weight"])
            TimetotalsStratum$weight  <- weight_from_Swan_schedule
            if (weight_from_Swan_schedule != 1) { cat("Note: weights are derived from Swan_schedule and != 1", "\n") }
          }
          
          #         add to total file
          TimetotalsAllStrata <-rbind(TimetotalsAllStrata, TimetotalsStratum) # combines TT-files by putting the time totals als counts in one file
          
        } # end loop strata
        write.table(TimetotalsAllStrata,paste0(working_folder,"\\",all_species[s], "_",  level-1, "_","TimetotalsAllStrata.csv"), dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
        
        colnames(TimetotalsAllStrata) <- c("time", "site", "count", "weight")  # RTRIM requires terms site and time
        
        ###########################################################################
        #       combine ocv files of all strata of pool together                  #
        ###########################################################################
        
        icv = list() # defines file in list format
        
        for (st in 1:nstrata) {
          
          if (level == 2) { 
            filenameCovmatrix <- paste(all_species[s], "_",  level-1, "_",  strata[st], "_ocv.csv", sep="")  # level1 = country file input
          } 
          if (level > 2) { 
            filenameCovmatrix <- paste(schedule_code, "_", all_species[s], "_", level-1, "_", strata[st], "_ocv.csv", sep="")   # pool file input, with schedule in name  
          }
          #print(paste0("Leyendo ",filenameCovmatrix))
          icv_1site = unname(as.matrix(read.table(filenameCovmatrix,header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\""))) # unname is required to delete colnames and rownames #13/02/2023 Javi: change "," for "." as decimal symbol 
          icv[[st]] = icv_1site # combine ocv's in a list needed for RTRIM
          
        } # end loop strata
        
        ###########################################################################
        #       run RTRIM using tt and icv of all strata combined                 #
        ###########################################################################
        
        weight                  <- TimetotalsAllStrata$weight
        
        result <- tryCatch(trim(count ~ site + time, data=TimetotalsAllStrata, weights = "weight", model = 2, changepoints = "all", autodelete = TRUE, covin=icv),error=warning) 
        #       default values for running RTRIM used
        
        if (class(result) == "trim") { 
          
          cat("Indices successfully produced for species", all_species[s], "Level", level, "Pool", pools[p], "\n")
          save(x = result, file = paste(schedule_code, "_", all_species[s], "_", level, "_", pools[p],".RData", sep = ""))  # save as RData
          
          #         tt output (in a format simlar to that produced by RTRIM shell) 
          tt <- totals(result, which = "both")            
          filenameTimeTotals_out <- paste(schedule_code, "_",all_species[s], "_",  level, "_",  pools[p], "_indices_TT.csv", sep="")   # save as output to be used in next leve 
          #print(paste0("Time: ",TimetotalsAllStrata$time," base_year =", base_year  ))
          if(min(TimetotalsAllStrata$time)< base_year){        #Javi 21.09.2021
            indices <- index(result, which = "both",base= base_year)
          }else{
            indices <- index(result, which = "both")
            }                                             #Javi 21.09.2021
          indices$time <- NULL                     
          tt_indices <- cbind(tt, indices)         
          
          #         RTRIM requires terms site and time, and colnames should be renamed again as Year en TT_imputed
          colnames(tt_indices) <- c("Year", "TT_model", "TT_model_SE", "TT_imputed", "TT_imputed_SE", "Index_model", "Index_model_SE", "Index_imputed", "Index_imputed_SE") 
          write.table(tt_indices, file=filenameTimeTotals_out, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
          
          # ARCO 9-10-2021
          # File with fitted values
          FI <- results(result)
          # storage
          filenameFittedValues_out <- paste(schedule_code, "_",all_species[s], "_",  level, "_",  pools[p], "_fitted_values.csv", sep="")   # save as output for information 
          write.table(FI, filenameFittedValues_out, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
          # ARCO 9-10-2021    
          
          
          #         ocv output that may be used in the next level  # similar to RTRIM shell output
          covmatrix <- vcov(result)
          filenameOcv_out <- paste(schedule_code,"_", all_species[s], "_",  level, "_",  pools[p], "_ocv.csv", sep="")   # save as output to be used in next level
          write.table(covmatrix,file=filenameOcv_out, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
          
          #         fill arg_output_file # this file makes the output similar to output of RTRIM shell 
          
          if (exists("arg_output_file") == FALSE)  { cat("Error: arg_out file is missing; include Make_empty_arg_output_file in function list", "\n") }
          
          if (exists("arg_output_file") == TRUE)  {
            
            arg_output_file$N_sites                      <- result$nsite               # Number of unique sites in the counts file used. 
            arg_output_file$N_time_values                <- result$ntime               # Number of unique years in the counts file used.
            counts                                       <- tt_indices[, "TT_imputed"]
            arg_output_file$N_observed_zero_counts       <- sum(counts == 0)           # Number of zero counts in the entire period, in the counts file used.
            arg_output_file$N_observed_positive_counts   <- sum(counts > 0)            # Number of positive counts in the entire period, in the counts file used.
            arg_output_file$N_missing_counts             <- sum(is.na(counts))         # Number of missing counts in the entire period, in the counts file used.
            arg_output_file$N_counts_total               <- sum(counts > 0)              # Total number of counts in the entire period, in the counts file used.
            arg_output_file$Base_year_first_year         <- tt_indices[tt_indices$Index_imputed == 1 & tt_indices$Index_imputed_SE == 0, "Year"]  # for case of accidental more occurences of 100% Arco 16-2-2021
            arg_output_file$Base_year_last_year          <- tt_indices[tt_indices$Index_imputed == 1 & tt_indices$Index_imputed_SE == 0, "Year"]  # for case of accidental more occurences of 100% Arco 16-2-2021
            arg_output_file$Changepoints                 <- paste(result$changepoints, collapse = "-")    # Changepoints used to calculate the indices. #29/03/2023 Javi Change points shall have - instead , as separator, function modified accordingly.
            
            slopes_imputed <- overall(result, which = "imputed")
            arg_output_file$Slope_imputed_mul            <- slopes_imputed$slope$mul                       # Multiplicative imputed slope over the entire period.
            arg_output_file$Slope_imputed_mul_SE         <- slopes_imputed$slope$se_mul                    # Standard error of multiplicative imputed slope over the entire period.
            arg_output_file$Slope_imputed_classification <- slopes_imputed$slope$meaning                   # Trend classification of multiplicative imputed slope over the entire period.
            # overdisperison, serial correlation and slopes_from were not computed; all fields remain value 0 although NA would be more appropriate 
            filenameArg_out <- paste(schedule_code,"_", all_species[s], "_",  level, "_",  pools[p], "_arg_output.csv", sep="")   # save as output to be used in next level
            write.table(arg_output_file, filenameArg_out, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
          } 
          
        } else {
          
          cat("Failed to produce indices for species", all_species[s], "Level", level, "Pool", pools[p], "\n")
        }
        
      } # end loop pools
      
    } # end loop species
    
  } # end nlevels
  
} # E N D of function Rtrimcovin
##########################################################################################################################################################








Assess_range_of_years_for_resultsfile <- function(no_arguments) {  
  ####################################################################################
  # Aim: assess min and max year for later use in results files                      #
  # Input: indices_TT.csv files of all active Level1 strata                          # 
  # Output: Message in Logfile                                                       # 
  ####################################################################################
  
  cat("\n")
  cat("@@@ function Assess_range_of_years_for_resultsfile", "\n") 
  cat("\n") 
  
  Swan_active_records <- Swan_schedule_active[Swan_schedule_active$B_sel == 1, ]
  nrecords            <- nrow(Swan_active_records)
  
  allyear <- NULL 
  
  for (st in 1:nrecords) { 
    
    species_nr <- Swan_active_records[st, "Species_nr"]
    stratum_nr <- Swan_active_records[st, "Level1"]
    
    filenameTimeTotals_out <- paste(species_nr, "_",  1, "_",  stratum_nr, "_indices_TT.csv", sep="")  # level1 = country files input; no schedule_code in the name 
    TimetotalsStratum  <- read.table(filenameTimeTotals_out, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
    stratum_year       <- unique(TimetotalsStratum[, "Year"])
    allyear            <- c(allyear, stratum_year) 
    
  } # end loop strata
  
  ###########################################################################
  # store lowest and highest year for later use in another function         #
  ###########################################################################
  
  lowest_year_of_schedule  <<- min(allyear)
  highest_year_of_schedule <<- max(allyear)
  cat("Lowest_year_of_schedule currently is : ", lowest_year_of_schedule, "\n") 
  cat("Highest_year_of_schedule currently is: ", highest_year_of_schedule, "\n") 
  cat("\n")
  
} # E N D of function Assess_range_of_years_for_resultsfile
##########################################################################################################################################################














Gather_results_of_countries_horizontal_format  <- function(no_arguments) {
  ############################################################################
  # Aim: Collect output RTRIM from countries selected in active schedule     #
  # Input: Indices_TT & arg_output files (not RData)                         #
  # Output: Swan-out with indices and trends of active schedule              #
  ############################################################################
  
  cat("\n") 
  cat("@@@ function Gather_results_of_countries_horizontal_format executed", "\n")
  cat("\n") 
  
  ###########################################################################
  # assess number of records in results file                                # 
  # = combinations of species * countries in active schedule                #
  ###########################################################################
  
  Swan_schedule_active$combination_column <-paste(as.character(Swan_schedule_active[,"Species_nr"]), as.character(Swan_schedule_active[,"Level1"]), sep="-") 
  total_number_records    <- length(unique(Swan_schedule_active[, "combination_column"]))
  
  uyear <- lowest_year_of_schedule:highest_year_of_schedule
  nyear <- length(uyear)
  
  Swan_ind <- array(NA,dim=c(total_number_records, (18+nyear)))  
  colnames(Swan_ind) <- c("Schedule_code", "Schedule_name", "Stratum_type", "Stratum_number", "Stratum_name", "Year_of_analysis", "Species_number", "Recordtype_number", "Recordtype_name", "N_sites", "Slope_imputed_mul", "Slope_imputed_mul_SE", "Slope_imputed_classification", "Slope_from_imputed_mul", "Slope_from_imputed_mul_SE", "Slope_from_imputed_classification", "Year_from", "Date_analysis", uyear) 
  Swan_ind = data.frame(Swan_ind, check.names = FALSE) # check.names = FALSE prevents setting X before calendar year as colnames
  Swan_ind_se    <- Swan_ind
  Swan_tt        <- Swan_ind
  Swan_tt_se     <- Swan_ind
  
  ###########################################################################
  # loops species * countries                                               #
  ###########################################################################
  # Note: loops of species and countries here instead of loop of records from Swan_schedule; else the same species & countrynames need to be handeld multiple times
  
  species_stratum_rec <- 1  # record counter
  
  all_species <- unique(Swan_schedule_active[Swan_schedule_active$B_sel == 1, "Species_nr"])
  nspecies    <- length(all_species)
  
  for (s in 1:nspecies) { 
    
    Species      <- all_species[s]
    Species_name <- as.character(species_file[species_file$Species_nr == Species, "Species_name"] )
    print(Species)
    countries    <- unique(Swan_schedule_active[Swan_schedule_active$Species_nr == all_species[s], "Level1"]) # level1 contains country code
    ncountries   <- length(countries)
    
    for (c in 1:ncountries) { 
      
      Country      <- countries[c]
      print(Country)
      Country_name <- as.character(countries_file[countries_file$Country_code == Country, "Country"] ) 
      
      RTRIM_timetotals   <- paste(Species, "_1_", Country, "_indices_TT.csv",  sep="") 
      TimetotalsStratum  <- read.table(RTRIM_timetotals, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
      RTRIM_arguments    <- paste(Species, "_1_", Country, "_arg_output.csv",  sep="") 
      TimeargStratum     <- read.table(RTRIM_arguments, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
      nyear_stratum      <- TimeargStratum[, "N_time_values"]
      nplot_stratum      <- TimeargStratum[, "N_sites"]
      slope_imputed_mul            <- TimeargStratum[, "Slope_imputed_mul"]
      slope_imputed_mul_se         <- TimeargStratum[, "Slope_imputed_mul_SE"]
      slope_imputed_classification <- TimeargStratum[, "Slope_imputed_classification"]
      
      #     Slope_from ignored as Year_from varies in the files which countries have provided
      
      ###########################################################################
      #     fill in tables                                                      #
      ###########################################################################
      
      #     fill in table with indices
      
      Swan_ind[species_stratum_rec, "Schedule_code"]              <- schedule_code
      Swan_ind[species_stratum_rec, "Schedule_name"]              <- schedule_name
      Swan_ind[species_stratum_rec, "Stratum_type"]               <- 1
      Swan_ind[species_stratum_rec, "Stratum_number"]             <- Country
      Swan_ind[species_stratum_rec, "Stratum_name"]               <- Country_name
      Swan_ind[species_stratum_rec, "Year_of_analysis"]           <- max(TimetotalsStratum[, "Year"])  # Year_of_analysis = most recent year with data
      Swan_ind[species_stratum_rec, "Species_number"]             <- Species
      Swan_ind[species_stratum_rec, "Recordtype_number"]          <- 1  
      Swan_ind[species_stratum_rec, "Recordtype_name"]            <- "indices"
      Swan_ind[species_stratum_rec, "N_sites"]                    <- nplot_stratum
      
      #     fill in table with indices se's
      Swan_ind_se[species_stratum_rec, "Schedule_code"]           <- schedule_code
      Swan_ind_se[species_stratum_rec, "Schedule_name"]           <- schedule_name
      Swan_ind_se[species_stratum_rec, "Stratum_type"]            <- 1
      Swan_ind_se[species_stratum_rec, "Stratum_number"]          <- Country
      Swan_ind_se[species_stratum_rec, "Stratum_name"]            <- Country_name
      Swan_ind_se[species_stratum_rec, "Year_of_analysis"]        <- max(TimetotalsStratum[, "Year"])
      Swan_ind_se[species_stratum_rec, "Species_number"]          <- Species
      Swan_ind_se[species_stratum_rec, "Recordtype_number"]       <- 2   
      Swan_ind_se[species_stratum_rec, "Recordtype_name"]         <- "se_indices"
      Swan_ind_se[species_stratum_rec, "N_sites"]                 <- nplot_stratum
      
      #     fill in table with time totals
      Swan_tt[species_stratum_rec, "Schedule_code"]               <- schedule_code
      Swan_tt[species_stratum_rec, "Schedule_name"]               <- schedule_name
      Swan_tt[species_stratum_rec, "Stratum_type"]                <- 1
      Swan_tt[species_stratum_rec, "Stratum_number"]              <- Country
      Swan_tt[species_stratum_rec, "Stratum_name"]                <- Country_name
      Swan_tt[species_stratum_rec, "Year_of_analysis"]            <- max(TimetotalsStratum[, "Year"])
      Swan_tt[species_stratum_rec, "Species_number"]              <- Species
      Swan_tt[species_stratum_rec, "Recordtype_number"]           <- 3
      Swan_tt[species_stratum_rec, "Recordtype_name"]             <- "time_totals"
      Swan_tt[species_stratum_rec, "N_sites"]                     <- nplot_stratum
      
      #     fill in table with time totals se's
      Swan_tt_se[species_stratum_rec, "Schedule_code"]            <- schedule_code
      Swan_tt_se[species_stratum_rec, "Schedule_name"]            <- schedule_name
      Swan_tt_se[species_stratum_rec, "Stratum_type"]             <- 1
      Swan_tt_se[species_stratum_rec, "Stratum_number"]           <- Country
      Swan_tt_se[species_stratum_rec, "Stratum_name"]             <- Country_name
      Swan_tt_se[species_stratum_rec, "Year_of_analysis"]         <- max(TimetotalsStratum[, "Year"])
      Swan_tt_se[species_stratum_rec, "Species_number"]           <- Species
      Swan_tt_se[species_stratum_rec, "Recordtype_number"]        <- 4
      Swan_tt_se[species_stratum_rec, "Recordtype_name"]          <- "se_time_totals"
      Swan_tt_se[species_stratum_rec, "N_sites"]                  <- nplot_stratum
      
      #     select the proper column of calendar years taking into account that countries have different years in their data
      for (j in 1:nyear_stratum) { 
        target_year <- TimetotalsStratum[j, "Year"]  
        col_year    <- which(colnames(Swan_ind) == target_year) # column is similar in all Swan output arrays
        Swan_ind[species_stratum_rec, col_year]    <- round(TimetotalsStratum[j, "Index_imputed"] * 100)
        Swan_ind_se[species_stratum_rec, col_year] <- round(TimetotalsStratum[j, "Index_imputed_SE"] * 100)
        Swan_tt[species_stratum_rec, col_year]     <- round(TimetotalsStratum[j, "TT_imputed"])
        Swan_tt_se[species_stratum_rec, col_year]  <- round(TimetotalsStratum[j, "TT_imputed_SE"])
      }  
      
      #     overall slope entire period; possibly only multiplicative slope available in country output 
      Swan_ind[species_stratum_rec, "Slope_imputed_mul"]            <- slope_imputed_mul
      Swan_ind[species_stratum_rec, "Slope_imputed_mul_SE"]         <- slope_imputed_mul_se
      Swan_ind[species_stratum_rec, "Slope_imputed_classification"] <- as.character(slope_imputed_classification)
      
      #     Slope_from ignored as Year_from varies between countries
      Swan_ind[species_stratum_rec, "Slope_from_imputed_mul"]             <- NA
      Swan_ind[species_stratum_rec, "Slope_from_imputed_mul_SE"]          <- NA
      Swan_ind[species_stratum_rec, "Slope_from_imputed_classification"]  <- NA
      Swan_ind[species_stratum_rec, "Year_from"]                          <- NA
      
      species_stratum_rec <- species_stratum_rec + 1 
      
    } # end loop countries
    
  } # end loop species
  
  cat("All species*country files processed", "\n")
  
  # write output in horizontal format in working_folder, not in output folder
  Swan_total                    <- rbind(Swan_ind, Swan_ind_se, Swan_tt, Swan_tt_se) 
  Swan_total[, "Date_analysis"] <- Sys.Date()
  Swan_total                    <- Swan_total[order(Swan_total[, "Species_number"], Swan_total[,"Stratum_number"], Swan_total[,"Recordtype_number"], decreasing=FALSE),] # Stratum_type always 1
  Swan_countries_name           <- paste(schedule_code, "_countries_all_Indices_All_Trends.csv", sep="")  
  
  write.table(Swan_total, file = Swan_countries_name, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  
  cat("all_Indices_All_trends.csv output produced for countries in schedule", schedule_name, "\n")
  
} # E N D of function Gather_results_of_countries_horizontal_format
##########################################################################################################################################################




Gather_results_of_countries_vertical_format  <- function(no_arguments) {
  ############################################################################
  # Aim: Collect output RTRIM from countries selected in active schedule     #
  # Input: Indices_TT & arg_output files (not RData)                         #
  # Output: Swan-out with indices and time totals                            #
  ############################################################################
  
  cat("\n") 
  cat("@@@ function Gather_results_of_countries_vertical_format executed", "\n")
  cat("\n") 
  
  ###########################################################################
  # loops species * countries                                               #
  ###########################################################################
  # Note: loops of species and countries here instead of loop of records from Swan_schedule; else the same species & countrynames need to be handled multiple times
  
  all_species <- unique(Swan_schedule_active[Swan_schedule_active$B_sel == 1, "Species_nr"])
  nspecies    <- length(all_species)
  
  trimind_countries <- NULL
  trimtot_countries <- NULL
  
  for (s in 1:nspecies) { 
   
    Species      <- all_species[s]
    Species_name <- as.character(species_file[species_file$Species_nr == Species, "Species_name"] )
    
    countries    <- unique(Swan_schedule_active[Swan_schedule_active$Species_nr == all_species[s], "Level1"]) # level1 contains country code
    ncountries   <- length(countries)
    
    for (c in 1:ncountries) { 
      
      Country      <- countries[c]
      
      Country_name <- as.character(countries_file[countries_file$Country_code == Country, "Country"] ) 
      
      RTRIM_timetotals   <- paste(Species, "_1_", Country, "_indices_TT.csv",  sep="") 
      TimetotalsStratum  <- read.table(RTRIM_timetotals, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
      RTRIM_arguments    <- paste(Species, "_1_", Country, "_arg_output.csv",  sep="") 
      TimeargStratum     <- read.table(RTRIM_arguments, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
      #     slope_imputed_mul            <- TimeargStratum[, "Slope_imputed_mul"]
      #     slope_imputed_mul_se         <- TimeargStratum[, "Slope_imputed_mul_SE"]
      #     slope_imputed_classification <- TimeargStratum[, "Slope_imputed_classification"]
      
      #     Slope_from ignored as Year_from may vary between countries
      
      ###########################################################################
      #     fill in tables                                                      #
      ###########################################################################
      
      TimetotalsStratum$Code          <- all_species[s]
      TimetotalsStratum$Species       <- Species_name
      TimetotalsStratum$StratumType   <- 1
      TimetotalsStratum$Stratum       <- Country
      TimetotalsStratum$CountryGroup  <- Country_name
      
      TimetotalsStratum$Index         <- TimetotalsStratum$Index_imputed
      TimetotalsStratum$Index_SE      <- TimetotalsStratum$Index_imputed_SE
      trimind_part                    <- subset(TimetotalsStratum, select = c(Code, Species, StratumType, Stratum, CountryGroup, Year, Index, Index_SE))
      trimind_countries               <- rbind(trimind_countries, trimind_part)
      
      TimetotalsStratum$TT_Model      <- TimetotalsStratum$TT_model
      TimetotalsStratum$TT_Model_SE   <- TimetotalsStratum$TT_model_SE
      TimetotalsStratum$TT_Imputed    <- TimetotalsStratum$TT_imputed
      TimetotalsStratum$TT_Imputed_SE <- TimetotalsStratum$TT_imputed_SE
      
      trimtot_part                    <- subset(TimetotalsStratum, select = c(Code, Species, StratumType, Stratum, CountryGroup, Year, TT_Model, TT_Model_SE, TT_Imputed, TT_Imputed_SE))
      trimtot_countries               <- rbind(trimtot_countries, trimtot_part) 
      
    } # end loop countries
    
  } # end loop species
  
  # delete any double records
  trimind_countries$Combicolumn <- paste(as.character(trimind_countries[, "Code"]), as.character(trimind_countries[, "Stratum"]), as.character(trimind_countries[, "Year"]), sep="-") 
  trimind_countries             <- subset(trimind_countries, !duplicated(Combicolumn))    
  trimind_countries$Combicolumn <- NULL
  
  trimtot_countries$Combicolumn <- paste(as.character(trimtot_countries[, "Code"]), as.character(trimtot_countries[, "Stratum"]), as.character(trimtot_countries[, "Year"]), sep="-") 
  trimtot_countries             <- subset(trimtot_countries, !duplicated(Combicolumn))    
  trimtot_countries$Combicolumn <- NULL
  
  cat("All species*country files processed", "\n")
  
  # write output in vertical format in working_folder, not in output folder
  write.table(trimind_countries, file = "Trimind_countries.csv",dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  write.table(trimtot_countries, file = "Trimtot_countries.csv", dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  
  cat("trimtot.csv and trimind.csv produced for countries in schedule in vertical format in schedule", schedule_name, "\n")
  
} # E N D of function Gather_results_of_countries_vertical_format
##########################################################################################################################################################









Gather_results_of_pools_horizontal_format <- function(years) { 
  ############################################################################
  # Aim: Collect output RTRIM from pools selected in active schedule         #
  # Input: RData (not Indices_TT & arg_output files)                         #
  # Arguments with base_year and slopes_from                                 #
  # Output: Swan-out with indices and trends of active schedule              #
  # Output: Swan-out with slopes_from of active schedule in another format   #
  # All pools have the same base year and slopes-from                        #
  ############################################################################
  
  cat("\n") 
  cat("@@@ function Gather_results_of_pools_horizontal_format executed", "\n")
  cat("\n") 
  
  ###########################################################################
  # read arguments                                                          #
  ###########################################################################
  
  base_year        <- as.numeric(years[1])
  cat("base_year is", base_year, "\n")
  
  years_slope_from <- NULL
  for (y in 2:length(years)) {
    if (!is.na(years[y])) { 
      years_slope_from <- c(years_slope_from, years[y])
    }
  } 
  cat("years with slopes_from are", years_slope_from, "\n")
  
  ###########################################################################
  # assess number of records in resuls file                                 # 
  # = combinations of level * species * pool  in Swan_schedules             #
  ###########################################################################
  total_number_records <- 0
  
  LevelCols <- as.vector(grep("Level", names(Swan_schedules)))
  
  # look up number of levels in Swan_schedules.csv (similar for all schedules)
  nlevels <- length(LevelCols)  
  
  for (level in 2:nlevels) { 
    #   combine cols with species number
    Swan_schedule_active$Combikolom <- paste(as.character(Swan_schedule_active$Species), as.character (Swan_schedule_active[, LevelCols[level]] ) ,sep="") 
    number_records                  <- length(unique(Swan_schedule_active$Combikolom))
    total_number_records            <- total_number_records + number_records
  }
  
  uyear <- lowest_year_of_schedule:highest_year_of_schedule # obtained in earlier function
  last_year_in_data <- max(uyear)
  #print(paste0("Este es el last year in data ",last_year_in_data))
  nyear <- length(uyear)
  
  Swan_ind <- array(NA,dim=c(total_number_records, (18+nyear)))  
  colnames(Swan_ind) <- c("Schedule_code", "Schedule_name", "Stratum_type", "Stratum_number", "Stratum_name", "Year_of_analysis", "Species_number", "Recordtype_number", "Recordtype_name", "N_sites", "Slope_imputed_mul", "Slope_imputed_mul_SE", "Slope_imputed_classification", "Slope_from_imputed_mul", "Slope_from_imputed_mul_SE", "Slope_from_imputed_classification", "Year_from", "Date_analysis", uyear) 
  Swan_ind = data.frame(Swan_ind, check.names = FALSE) # check.names = FALSE prevents setting X before calendar year as colnames
  Swan_ind_se    <- Swan_ind
  Swan_tt        <- Swan_ind
  Swan_tt_se     <- Swan_ind
  
  # separate Swan with all slopes_from (i.e. for each Swan_ind record several slopes from may be calculated)
  total_number_records_slopes <- total_number_records * length(years_slope_from)
  trimslo_from <- array(NA,dim=c(total_number_records_slopes, 11)) 
  colnames(trimslo_from)    <- c("Code", "Species", "StratumType", "Stratum", "CountryGroup", "Slope Add", "Slope Add SE", "Slope Mult", "Slope Mult SE", "Trend Classification", "Year from")
  trimslo_from = data.frame(trimslo_from, check.names = FALSE) # check.names = FALSE prevents setting X before calendar year as colnames
  
  ###########################################################################
  # loops levels * species * pools                                          #
  ###########################################################################
  
  species_stratum_rec        <- 1  # record counter for Swan-output with indices
  species_stratum_slopes_rec <- 1  # record counter for Swan-output with slopes_from
  
  all_species <- unique(Swan_schedule_active[Swan_schedule_active$B_sel == 1, "Species_nr"])
  nspecies    <- length(all_species)
  
  for (level in 2:nlevels) { # loop over levels 
    
    Pool_column                     <- paste("Level", level, sep ="") 
    #print(Pool_column)
    Pool_column_before              <- paste("Level", level - 1, sep ="") 
    
    for (s in 1:nspecies) { 
      
      Species_name <- as.character(species_file[species_file$Species_nr == all_species[s],"Species_name"])
      #print(Species_name)
      pools     <- unique(Swan_schedule_active[Swan_schedule_active$Species == all_species[s], Pool_column])
      #print(pools)
      npools    <- length(pools)
      
      for (p in 1:npools) { 
        
        Pool_name <- as.character(pools_file[pools_file$Pool == pools[p], "Pool_Description"])
        #print(pools_file$Pool == pools[p])
        file <- paste(schedule_code, "_", all_species[s], "_", level, "_", pools[p],".RData", sep = "") 
        #print(file)
        if (file.exists(file)  == FALSE) {cat("Error: missing file", file,"\n")   }  
        
        if (file.exists(file)  == TRUE) { 
          
          load(file) # rtrim result with name as saved 
          
          #         retrieve information from RData
          #print("doing the model")
          tt                 <- totals(result, which = c("imputed")) 
          first_year         <- tt[1, "time"] 
          if (base_year < first_year)  { base_year_ok <- first_year  }
          if (base_year >= first_year) { base_year_ok <- base_year   }
          indices            <- index(result, base = base_year_ok) 
          nyear_stratum      <- nrow(tt)
          #print("we did the model")
          #         filter slopes_from that match with the years available in the data                                     # Arco 19-12-2020
          #print("Filtering the slopes")
          years_in_data     <- tt[, "time"]                                                                                # Arco 19-12-2020
          years_in_data     <- years_in_data[1:nrow(tt)-1]  # slopes_from not possible for last year, so skip last year    # Arco 19-12-2020
          years_filtered    <- is.element(years_slope_from, years_in_data)                                                 # Arco 19-12-2020
          years_slope_from_filtered <- years_slope_from[years_filtered]                                                    # Arco 19-12-2020
          #print("slopes filtered")
          ###########################################################################
          #         fill in tables                                                  #
          ###########################################################################
          
          #         fill in table with indices se's
          #print(Swan_ind[species_stratum_rec,c(1:15)])
          Swan_ind[species_stratum_rec, "Schedule_code"]          <- schedule_code
          #print(paste0("zzz",schedule_code))
          Swan_ind[species_stratum_rec, "Schedule_name"]          <- schedule_name
          #print(paste0("xxx",schedule_name))
          Swan_ind[species_stratum_rec, "Stratum_type"]           <- level
          #print(paste0("ccc",level))
          Swan_ind[species_stratum_rec, "Stratum_number"]         <- pools[p]
          #print(paste0("vvv",pools[p]))
          Swan_ind[species_stratum_rec, "Stratum_name"]           <- Pool_name
          #print(paste0("bbb",Pool_name))
          Swan_ind[species_stratum_rec, "Year_of_analysis"]       <- last_year_in_data
          #print(paste0("nnn",last_year_in_data))
          Swan_ind[species_stratum_rec, "Species_number"]         <- all_species[s]
          #print(all_species[s])
          Swan_ind[species_stratum_rec, "Recordtype_number"]      <- 1
          Swan_ind[species_stratum_rec, "Recordtype_name"]        <- "indices"
          Swan_ind[species_stratum_rec, "N_sites"]                <- result$nsite # number of strata in the pool
          #print("we filled the indeces")
          #         fill in table with indices se's
          Swan_ind_se[species_stratum_rec, "Schedule_code"]       <- schedule_code
          Swan_ind_se[species_stratum_rec, "Schedule_name"]       <- schedule_name
          Swan_ind_se[species_stratum_rec, "Stratum_type"]        <- level
          Swan_ind_se[species_stratum_rec, "Stratum_number"]      <- pools[p]
          Swan_ind_se[species_stratum_rec, "Stratum_name"]        <- Pool_name
          Swan_ind_se[species_stratum_rec, "Year_of_analysis"]    <- last_year_in_data
          Swan_ind_se[species_stratum_rec, "Species_number"]      <- all_species[s]
          Swan_ind_se[species_stratum_rec, "Recordtype_number"]   <- 2  
          Swan_ind_se[species_stratum_rec, "Recordtype_name"]     <- "se_indices"
          Swan_ind_se[species_stratum_rec, "N_sites"]             <- result$nsite # number of strata in the pool
          #print("Wi filled indeces se's")
          #         fill in table with time totals 
          Swan_tt[species_stratum_rec, "Schedule_code"]           <- schedule_code
          Swan_tt[species_stratum_rec, "Schedule_name"]           <- schedule_name
          Swan_tt[species_stratum_rec, "Stratum_type"]            <- level
          Swan_tt[species_stratum_rec, "Stratum_number"]          <- pools[p]
          Swan_tt[species_stratum_rec, "Stratum_name"]            <- Pool_name
          Swan_tt[species_stratum_rec, "Year_of_analysis"]        <- last_year_in_data
          Swan_tt[species_stratum_rec, "Species_number"]          <- all_species[s]
          Swan_tt[species_stratum_rec, "Recordtype_number"]       <- 3
          Swan_tt[species_stratum_rec, "Recordtype_name"]         <- "time_totals"
          Swan_tt[species_stratum_rec, "N_sites"]                 <- result$nsite # number of strata in the pool
          #print("We filled the TT's")
          #         fill in table with time totals se's
          Swan_tt_se[species_stratum_rec, "Schedule_code"]        <- schedule_code
          Swan_tt_se[species_stratum_rec, "Schedule_name"]        <- schedule_name
          Swan_tt_se[species_stratum_rec, "Stratum_type"]         <- level
          Swan_tt_se[species_stratum_rec, "Stratum_number"]       <- pools[p]
          Swan_tt_se[species_stratum_rec, "Stratum_name"]         <- Pool_name
          Swan_tt_se[species_stratum_rec, "Year_of_analysis"]     <- last_year_in_data
          Swan_tt_se[species_stratum_rec, "Species_number"]       <- all_species[s]
          Swan_tt_se[species_stratum_rec, "Recordtype_number"]    <- 4
          Swan_tt_se[species_stratum_rec, "Recordtype_name"]      <- "se_time_totals"
          Swan_tt_se[species_stratum_rec, "N_sites"]              <- result$nsite # number of strata in the pool
          #print("We filled TT's SE")
          
          #         select the proper column of calendar years taking into account that countries have different years in their data
          for (j in 1:nyear_stratum) { 
            target_year <- tt[j, "time"]  
            col_year    <- which(colnames(Swan_ind) == target_year) # column is similar in alll Swan arrays
            Swan_ind[species_stratum_rec, col_year]               <- round(indices[j, "imputed"] * 100)
            Swan_ind_se[species_stratum_rec, col_year]            <- round(indices[j, "se_imp"] * 100)
            Swan_tt[species_stratum_rec, col_year]                <- tt[j, "imputed"]
            Swan_tt_se[species_stratum_rec, col_year]             <- tt[j, "se_imp"]
          }  
          
          #         overall multiplicative slope of entire periode 
          #print("Gathering the slopes")
          slope <- overall(result)
          Swan_ind[species_stratum_rec, "Slope_imputed_mul"]          <- slope$slope$mul
          Swan_ind[species_stratum_rec, "Slope_imputed_mul_SE"]       <- slope$slope$se_mul
          Swan_ind[species_stratum_rec, "Slope_imputed_classification"] <- slope$slope$meaning  
          #print("Wegather the slopes_mul")
          ##########################################################################               # Arco 16-2-2021
          #        5820 Black-headed gull - a problematic species
          # AvS: Stripping of ocv files did not help me neither, so I had a closer look at what is causing the troubles. 
          #It occurs that RTRIM is unable to assess estimates of slopes_from using the data. 
          #The only solution I can think of is to skip calculating slopes_from for such species. 
          #It is not really a nice  solution, but at least you get indices and overall indices (see attached).
          ########################################################################## 
          
          # In order to achieve this, you have to skip part of the code by inserting the yellow lines in the function Gather_results_of_pools_horizontal_format.
          
          #if (all_species[s] != 5820) {                                                           # Arco 16-2-2021 - 1 line added
            #         slope_from year mentioned as first argument in function list # = must be: mentioned as second argument = first slope_from # Arco 1-7-2021
            
          possible_years_slope_from <- is.element(years_slope_from_filtered, tt[,1])                                             # Arco 27-11-2021
            #years_slope_from          <- years_slope_from[years_slope_from_possible]                                                  # Arco 1-7-2021
            
             #         slope_from year mentioned as first argument in function list
            if (possible_years_slope_from[1] == TRUE) {                                                                               # Arco 1-7-2021
            slope_from <- overall(result,changepoints=c(years_slope_from_filtered[1]))                                                            # Arco 19-12-2020
            Swan_ind[species_stratum_rec, "Slope_from_imputed_mul"]            <- slope_from$slope$mul
            Swan_ind[species_stratum_rec, "Slope_from_imputed_mul_SE"]         <- slope_from$slope$se_mul
            Swan_ind[species_stratum_rec, "Slope_from_imputed_classification"] <- slope_from$slope$meaning  
            Swan_ind[species_stratum_rec, "Year_from"]                         <- years_slope_from_filtered[1]                                    # Arco 19-12-2020 
            }                                                                                                                         # Arco 1-7-2021 
            #         fill in slopes_from in a separate outputfile with vertical format
            
            for (y in 1:length(years_slope_from_filtered)) {                                                                                      # Arco 19-12-2020   
              
              if (possible_years_slope_from[y] == TRUE) {                                                                             # Arco 1-7-2021
                trimslo_from[species_stratum_slopes_rec, "Code"]                   <- all_species[s]
                trimslo_from[species_stratum_slopes_rec, "Species"]                <- Species_name
                trimslo_from[species_stratum_slopes_rec, "StratumType"]            <- level
                trimslo_from[species_stratum_slopes_rec, "Stratum"]                <- pools[p]
                trimslo_from[species_stratum_slopes_rec, "CountryGroup"]           <- Pool_name
              
                slope_from <- overall(result,changepoints=c(years_slope_from_filtered[y]))                                                          # Arco 19-12-2020 
                trimslo_from[species_stratum_slopes_rec, "Slope Add"]              <- slope_from$slope$add
                trimslo_from[species_stratum_slopes_rec, "Slope Add SE"]           <- slope_from$slope$se_add
                trimslo_from[species_stratum_slopes_rec, "Slope Mult"]             <- slope_from$slope$mul
                trimslo_from[species_stratum_slopes_rec, "Slope Mult SE"]          <- slope_from$slope$se_mul
                trimslo_from[species_stratum_slopes_rec, "Trend Classification"]   <- slope_from$slope$meaning  
                trimslo_from[species_stratum_slopes_rec, "Year from"]              <- years_slope_from_filtered[y]                                  # Arco 19-12-2020   
              
                species_stratum_slopes_rec <- species_stratum_slopes_rec + 1  # counter for Swan-output with slopes_from
                #print("We gathered the slopes_mul")
              }                                                                                                                       # Arco 1-7-2021
            }
            
          #} # end of != 5820                                                             # Arco 16-2-2021 - 1 line added
          
          species_stratum_rec        <- species_stratum_rec + 1 # counter for Swan-output with indices
          
        } # if file exists          
        
      } # end loop pools
      
    } # end loop species
    
  } # end loop levels
  
  cat("All species*pools processed", "\n")
  
  # combine separate arrays  
  Swan_total <- rbind(Swan_ind, Swan_ind_se, Swan_tt, Swan_tt_se) 
  Swan_total[, "Date_analysis"]  <- Sys.Date()
  
  # combine files of pools and countries
  Swan_countries_name <-  paste(schedule_code, "_countries_all_Indices_All_Trends.csv", sep="")  # table with indices and se's for all species
  if (file.exists(Swan_countries_name))  { 
    Swan_countries <- read.table(Swan_countries_name, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
    Swan_total  <- rbind(Swan_total, Swan_countries) 
  } 
  
  Swan_total <- Swan_total[order(Swan_total[, "Species_number"], Swan_total[, "Stratum_type"], Swan_total[, "Stratum_number"], Swan_total[,"Recordtype_number"], decreasing=FALSE),]
  
  setwd(output_folder) 
  Swan_name <-  paste(schedule_code, "_all_Indices_All_Trends.csv", sep="")  # table with indices and se's for all species
  write.table(Swan_total, file = Swan_name, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol 
  cat("all_Indices_all_trends.csv output produced for pools and countries together in schedule", schedule_name, "\n")  
  
  #######################################################
  # write output slopes_from POOLS ONLY
  trimslo_from <- trimslo_from[!is.na(trimslo_from[,1]),]  # remove any empty record caused by slope_from years that were not possible     # Arco 1-7-2021
  trimslo_from <- trimslo_from[order(trimslo_from[, "Species"], trimslo_from[, "StratumType"], trimslo_from[,"Stratum"], decreasing=FALSE),]
  trimslo_name <- paste(schedule_code, "_", "Trimslo.csv", sep="")  # only pools. If countries had to be added, they first need to be converted into pools by running RSWAN for countries
  write.table(trimslo_from, file = trimslo_name, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  setwd(working_folder) 
  
  cat("trimslo.csv output produced for pools in schedule", schedule_name, "\n")
  
} # E N D of function Gather_results_of_pools_horizontal_format
##########################################################################################################################################################






Gather_results_of_pools_vertical_format <- function(years) {                              # Arco 19-12-2020  
  ############################################################################
  # Aim: Collect output RTRIM from pools selected in active schedule         #
  # Input: RData (not Indices_TT & arg_output files)                         #
  # Arguments with base_year and slopes_from                                 #
  # Output: Swan-out with indices and time totals                            #
  ############################################################################
  
  cat("\n") 
  cat("@@@ function Gather_results_of_pools_vertical_format executed", "\n")
  cat("\n") 
  
  ###########################################################################
  # read arguments                                                          #
  ###########################################################################
  
  base_year        <- as.numeric(years[1])
  cat("base_year is", base_year, "\n")
  
  ###########################################################################
  # loops levels * species * pools                                          #
  ###########################################################################
  
  LevelCols <- as.vector(grep("Level", names(Swan_schedules)))
  # look up number of levels in Swan_schedules.csv (similar for all schedules)
  nlevels <- length(LevelCols)  
  
  all_species <- unique(Swan_schedule_active[Swan_schedule_active$B_sel == 1, "Species_nr"])
  nspecies    <- length(all_species)
  
  trimind_pools <- NULL
  trimtot_pools <- NULL
  trimind_total <- NULL
  trimtot_total <- NULL
  
  for (level in 2:nlevels) { # loop over levels 
    
    Pool_column                     <- paste("Level", level, sep ="") 
    Pool_column_before              <- paste("Level", level - 1, sep ="") 
    
    for (s in 1:nspecies) { 
      
      Species_name <- as.character(species_file[species_file$Species_nr == all_species[s],"Species_name"])
      
      pools     <- unique(Swan_schedule_active[Swan_schedule_active$Species == all_species[s], Pool_column])
      npools    <- length(pools)
      
      for (p in 1:npools) { 
        
        Pool_name <- as.character(pools_file[pools_file$Pool == pools[p], "Pool_Description"])
        
        file <- paste(schedule_code, "_", all_species[s], "_", level, "_", pools[p],".RData", sep = "")
        if (file.exists(file)  == FALSE) {cat("Error: missing file", file,"\n")   }  
        
        if (file.exists(file)  == TRUE) { 
          
          load(file) # rtrim result with name as saved 
          
          #         retrieve information from RData
          tt                 <- totals(result, which = c("both")) 
          colnames(tt)       <- c("Year", "TT_Model", "TT_Model_SE", "TT_Imputed", "TT_Imputed_SE")
          first_year         <- tt[1, "Year"]     
          if (base_year < first_year)  { base_year_ok <- first_year  }
          if (base_year >= first_year) { base_year_ok <- base_year   }
          indices            <- index(result, base = base_year_ok) 
          colnames(indices)  <- c("Year", "Index", "Index_SE")
          nyear_stratum      <- nrow(tt)
          
          ###########################################################################
          #         fill in tables                                                  #
          ###########################################################################
          
          indices$Code         <- all_species[s]
          indices$Species      <- Species_name
          indices$StratumType  <- level
          indices$Stratum      <- pools[p]
          indices$CountryGroup <- Pool_name
          trimind_part         <- subset(indices, select = c(Code, Species, StratumType, Stratum, CountryGroup, Year, Index, Index_SE))
          trimind_pools        <- rbind(trimind_pools, trimind_part)
          
          tt$Code              <- all_species[s]
          tt$Species           <- Species_name
          tt$StratumType       <- level
          tt$Stratum           <- pools[p]
          tt$CountryGroup      <- Pool_name
          trimtot_part         <- subset(tt, select = c(Code, Species, StratumType, Stratum, CountryGroup, Year, TT_Model, TT_Model_SE, TT_Imputed, TT_Imputed_SE))
          trimtot_pools        <- rbind(trimtot_pools, trimtot_part) 
          
        } # if file exists          
        
      } # end loop pools
      
    } # end loop species
    
  } # end loop levels
  
  # delete any double records
  trimind_pools$Combicolumn <- paste(as.character(trimind_pools[, "Code"]), as.character(trimind_pools[, "Stratum"]), as.character(trimind_pools[, "Year"]), sep="-") 
  trimind_pools             <- subset(trimind_pools, !duplicated(Combicolumn))    
  trimind_pools$Combicolumn <- NULL
  
  trimtot_pools$Combicolumn <- paste(as.character(trimtot_pools[, "Code"]), as.character(trimtot_pools[, "Stratum"]), as.character(trimtot_pools[, "Year"]), sep="-") 
  trimtot_pools             <- subset(trimtot_pools, !duplicated(Combicolumn))    
  trimtot_pools$Combicolumn <- NULL
  
  cat("All species*pools processed", "\n")
  
  # combine output of pools and countries in one file
  
  if (file.exists("Trimind_countries.csv"))  { 
    trimind_countries <- read.table("Trimind_countries.csv", header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol  
    trimind_total     <- rbind(trimind_pools, trimind_countries)
  } 
  trimind_total <- trimind_total[order(trimind_total[, "Code"], trimind_total[, "StratumType"], trimind_total[,"Stratum"], decreasing=FALSE),]
  
  if (file.exists("Trimtot_countries.csv"))  { 
    trimtot_countries <- read.table("Trimtot_countries.csv",header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol  
    trimtot_total     <- rbind(trimtot_pools, trimtot_countries)
  } 
  trimtot_total <- trimtot_total[order(trimtot_total[, "Code"], trimtot_total[, "StratumType"], trimtot_total[,"Stratum"], decreasing=FALSE),]
  
  # write output in vertical format
  setwd(output_folder) 
  trimind_name <-  paste(schedule_code, "_Trimind.csv", sep="")  
  trimtot_name <-  paste(schedule_code, "_Trimtot.csv", sep="")  
  write.table(trimind_total, file = trimind_name,dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  write.table(trimtot_total, file = trimtot_name, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  setwd(working_folder) 
  
  cat("trimtot and trimind produced for pools and countries combined in vertical format", schedule_name, "\n")
  
} # E N D of function Gather_results_of_pools_vertical_format
##########################################################################################################################################################

Backup_files<- function (general_folder,working_folder){ #16/02/2023 Javi:Adding a new function to creat backup for general and working folders
  ############################################################################
  # Aim: Generate a backup folder for general and working folders            #
  # Input: general and working folder paths                                  #
  # Output: Backup_general and backup_working folders                        #
  ############################################################################
  cat("\n")
  cat("@@@ function backup_files executed", "\n")
  cat("\n")
  
  backup_general <-sub("general_folder\\\\", "backup_general_folder", general_folder)
  backup_working <-sub("working_folder\\\\", "backup_working_folder", working_folder)
  dir.create(backup_general)
  dir.create(backup_working)
  files_general=list.files(general_folder)
  files_working=list.files(working_folder)
  file.copy(from=paste0(general_folder,"\\",files_general),to=paste0(backup_general,"\\",files_general))
  file.copy(from=paste0(working_folder,"\\",files_working),to=paste0(backup_working,"\\",files_working))
  cat("\n")
  cat("Backup has been created", "\n")
  cat("\n")
}



Reset_rswan<- function (general_folder,working_folder){ #16/02/2023 Javi:Adding a new function to reset rswan to its original state ready to run again
  ############################################################################
  # Aim: Reset Rswan to use it again                                         #
  # Input: general and working folder paths                                  #
  # Output: repalces csv from the general folder for the general backup ones #
  # delete the working folder, creates it again and copy the files from the  #
  # working backup folder                                                    #
  ############################################################################
  cat("\n")
  cat(" Reseting Rswan", "\n")
  cat("\n")
  
  backup_general <-sub("general_folder\\\\", "backup_general_folder", general_folder)
  backup_working <-sub("working_folder\\\\", "backup_working_folder", working_folder)
  

  
  files_backup_general=list.files(backup_general,pattern = ".csv")
  file.copy(from=paste0(backup_general,"\\",files_backup_general),to=paste0(general_folder,"\\",files_backup_general),overwrite=TRUE)
  
  files_working=list.files(working_folder)
  files_backup_working=list.files(backup_working)
  
  unlink(working_folder,recursive = T)
  if(dir.exists(working_folder)==F){  
    dir.create(working_folder)}
  file.copy(from=paste0(backup_working,"\\",files_backup_working),to=paste0(working_folder,"\\",files_backup_working),overwrite=TRUE)
  cat("\n")
  cat("RSwan was reseted correctly", "\n")
  cat("\n")
}
#####################################################################################################################
#This function goes iterates through _indices_TT.csv _arg_outputs.csv and _ocv.csv files checking if their formating
#is correct. Focusing on the changes introduced to RTRIM where the decimal symbol was change to "." and the separator
#Of change points was changed to "-"
#####################################################################################################################
Check_input<-function(working_folder){ ##20/03/2023 Javi:Adding a new function to check the inputs formating
  #Checking arg_outputs files
  cat("\n")
  cat("Checking arg_ouput.csv files formating", "\n")
  cat("\n")
  arg_ouput_list<-list.files(working_folder,pattern = "_arg_output.csv")
  control=0
  for(i in 1:length(arg_ouput_list)){
    
    xx<- read.table(paste0(working_folder,arg_ouput_list[i]), header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
  
    if(grepl("-",xx$Changepoints,fixed=T)==F){
      
      cat(paste0("Changepoints separator in ",arg_ouput_list[i])," are different from '-' ", "\n")
      
      control=control+1
    }
    if(is.numeric(xx$Slope_imputed_mul)==F){
      
      cat(paste0(arg_ouput_list[i])," numerics arguments seem to have a different decymal symbol than '.'", "\n")
      
    }
  }
  
  if(control == 0){
    cat("\n")
    cat("Arg_output.csv files seems to have correct formating"," \n")
    cat("\n")
  }
  #Checking _indices_TT.csv files
  cat("\n")
  cat("Checking indices_TT.csv files formating", "\n")
  cat("\n")
  indices_TT_list<-list.files(working_folder,pattern = "_indices_TT.csv")
  control_tt=0
  
  for(i in 1:length(indices_TT_list)){
    
    xx<- read.table(paste0(working_folder,indices_TT_list[i]), header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
    
    if(is.character(xx$Index_model)==T){
      control_tt=control_tt+1
      
      cat(paste0(indices_TT_list[i])," numerics arguments seem to have a different decimal symbol than '.'", "\n")
      
    }
  }
  if(control_tt==0){
    cat("\n")
    cat("Indices_TT.csv files seems to have correct formating"," \n")
    cat("\n")
  }
  
  #Checking _ocv.csv files
  cat("\n")
  cat("Checking _ocv.csv files formating", "\n")
  cat("\n")
  ocv_list<-list.files(working_folder,pattern = "_ocv.csv")
  control_ocv=0
  
  for(i in 1:length(ocv_list)){
    
    xx<- read.table(paste0(working_folder,ocv_list[i]), header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
    
    if(is.character(xx[1])==T){
      control_ocv=control_tt+1
   
      cat(paste0(ocv_list[i])," numerics arguments seem to have a different decimal symbol than '.'", "\n")

    }
    if(nrow(xx)!=ncol(xx)){
      cat("\n")
      cat(paste0(ocv_list[i])," have wrong number of dimensions (ncol != nrow)", "\n")
      cat("\n") 
      control_ocv=control_ocv+1
    }
  }
  if(control_ocv==0){
    cat("\n")
    cat("ocv_list files seems to have correct formating"," \n")
    cat("\n")
  }
}
#END of the Inputs cheking function
  