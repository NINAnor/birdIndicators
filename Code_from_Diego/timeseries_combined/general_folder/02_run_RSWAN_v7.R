########################################################################################################################################
#Version V6 (script 02_run_RSWAN_V6=v5, but Pools and Schedules are changed according to changes in the database)

general_folder <- paste0(wd,"general_folder/")
working_folder <- paste0(wd,"working_folder/")
output_folder  <- paste(working_folder, "output/", sep = '') 
output_folder2 <- paste(working_folder, "output2/", sep = '') 

if (file.exists(output_folder2) == TRUE)  { unlink(output_folder2, recursive=TRUE) } 
if (file.exists(output_folder2) == FALSE) { dir.create(output_folder2) }  

setwd(general_folder)  

if (file.exists("03_Swan_Euromonitoring_v7.R") == FALSE) {cat("Error! R-code functions not found in general folder","\n")   }  

source ("03_Swan_Euromonitoring_v7.R")

Check_steering_files(general_folder)  

Run_Swan (general_folder, working_folder)  

###########################################################################################################################################
