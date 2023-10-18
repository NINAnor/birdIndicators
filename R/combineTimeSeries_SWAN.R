
combineTimeSeries_SWAN <- function(general_folder_abs, working_folder_abs){
  
  ## Ensure we reset to the project directory after executing this function
  project_wd <- getwd()
  on.exit(setwd(project_wd))
  
  ## Check whether required helper functions are available
  RSWAN_funNames <- c("Check_steering_files",
                      "Check_steering_files_combine",
                      "Run_Swan",
                      "Select_all_available_datafiles",
                      "Check_existence_datafiles",
                      "Check_content_datafiles",
                      "Check_missing_years_in_level2_pools",
                      "Import_weights_into_Swan_schedules",
                      "Truncate_datafiles",
                      "Strip_ocv_files",
                      "Make_empty_arg_output_file",
                      "Rtrimcovin",
                      "Assess_range_of_years_for_resultsfile",
                      "Gather_results_of_countries_horizontal_format",
                      "Gather_results_of_countries_vertical_format",
                      "Gather_results_of_pools_horizontal_format",
                      "Gather_results_of_pools_vertical_format",
                      "Backup_files",
                      "Reset_rswan"
                      )
  
  RSWAN_funExist <- rep(NA, length(RSWAN_funNames))
  
  for(i in 1:length(RSWAN_funNames)){
    RSWAN_funExist[i] <- exists(RSWAN_funNames[i])
  }
  
  if(any(!RSWAN_funExist)){
    missingFun <- RSWAN_funNames[which(!RSWAN_funExist)]
    stop(paste0("The following RSWAN helper functions are missing: ", 
                paste(missingFun, collapse = ", "),
                ". \nThese functions are provided in the script 03_Swan_Euromonitoring_v8.R"))
  }
  
  ## Add "/" to folder names for it to work with code
  general_folder <- paste0(general_folder_abs, "/")
  working_folder <- paste0(working_folder_abs, "/")
  
  # ----------------------------------------------------------------------------
  # START OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  ########################################################################################################################################
  #Version V6 (script 02_run_RSWAN_V6=v5, but Pools and Schedules are changed according to changes in the database)
  
  output_folder  <- paste(working_folder, "output/", sep = '') 
  output_folder2 <- paste(working_folder, "output", sep = '') 
  
  if (file.exists(output_folder2) == TRUE)  { unlink(output_folder2, recursive=TRUE) } 
  if (file.exists(output_folder2) == FALSE) { dir.create(output_folder2) }  
  
  #setwd(general_folder)  
  
  #if (file.exists("03_Swan_Euromonitoring_v8.R") == FALSE) {cat("Error! R-code functions not found in general folder","\n")   }  
  
  #source ("03_Swan_Euromonitoring_v8.R")
  ###########################################################################################################################################
  # Checking steering files and working folder files.
  
  # Checks steering files for the normal RSwan run
  Check_steering_files(general_folder) 
  
  # Checks steering files for the combine RSwan run
  Check_steering_files_combine(general_folder) 
  
  # Check inputs ensures that the formatting is correct, if formatting differs from the one given to the last RTRIM version
  # Files and a problem description will be printed in the console it checks all .csv files in the working folder.
  Check_input(working_folder)#Javi 20/03/2023
  
  ###########################################################################################################################################
  #Normal RSwan run
  
  
  Run_Swan (general_folder, working_folder)
  ###########################################################################################################################################
  #Combine the data files from old and new schemes (or different regions) into the final COMB_XX_XX_XX.csv files
  #Important: It requires a Swan_schedules_combine.csv telling which country (countricode-species needs to be combined
  # It also requires Species_Countries_combine.csv (copy of the original Species_Countries.csv) for the species-country population
  Run_Swan (general_folder, working_folder,1)
  ###########################################################################################################################################
  #Creates a backup for the general and working folder archives
  # IMPORTANT: Backup must be created before running RSWAN
  #Backup_files(general_folder,working_folder)
  ###########################################################################################################################################
  #Resets Rswan to its original state, importing the files from the backup folders
  # IMPORTANT: If there was some error in the Rswan execution, logfile might stop the function from working, close and open rswan before running
  # the reset function
  #Reset_rswan(general_folder,working_folder)

  # ----------------------------------------------------------------------------
  # END OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  message("The new RSWAN functions for making a backup and resetting RSWAN are not currently executed as they do not work if the names of the general_folder and working_folder are not exactly that.")
}
