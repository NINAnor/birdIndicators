#' Adds truncation of first survey year for specified species
#' 
#' This is a function wrapper for the PECBMS code file "01_CODE_adding_FirstSurvey_YEAR_OK.R"
#'
#' @param general_folder character. Path to the folder containing files 
#' provided by PECBMS for RSWAN analyses. This function requires the files
#' "Species_Countries.csv" and "New_notorics_2022.xlsx". 
#' @param working_folder character. Path to the working folder which contains
#' the PECBMS Trim results files for data from Hekkefuglovervåkingen and the 
#' predecessor monitoring.
#' @param maxYear integer. Latest year to include in analyses. 
#'
#' @return
#' @export
#'
#' @examples

correctFirstSurveyYear_SWAN <- function(general_folder, working_folder, maxYear){
  
  
  ## Add "/" to folder names for it to work with code
  general_folder <- paste0(general_folder, "/")
  working_folder <- paste0(working_folder, "/")
  
  # ----------------------------------------------------------------------------
  # START OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  species_country <- read.table(paste0(general_folder,"Species_Countries.csv"), header = T, check.names = FALSE, dec = ".", sep = ";", comment.char = "", quote = "\"")
  
  #Now we will add the exemption included in the table, for those species which normally needs to be truncated.
  #1st loading the data, be carefull with the route of the file.
  exemptions <- readxl::read_excel(paste0(general_folder,"New_notorics_2022.xlsx"))
  
  
  #LOOP: This loop currently uses the combination of species code and country code wich are in the same lane of the species_country file, this make the loop
  #more robust since is nor double nor dependent of the position of the data. Besides is simplier, and simplier is always good.
  
  
  for (i in 1:nrow(species_country)){
    
    try({
      x <- read.csv(paste0(working_folder, species_country$Species_nr[i], "_1_", species_country$Country_code[i], "_indices_TT.csv"), header = T, sep = ";")
      species_country$Year_first[i] <- min(x$Year)
      if(max(x$Year)>maxYear){
        species_country$Year_last[i] = maxYear
      }else{
        species_country$Year_last[i] <- max(x$Year)
      }
    }, silent = T)
  }
  
  
  #This double loop, look through both files, the table with exemptions and the species_country file, searching for a double match. First it matches the species code.
  #Second, it matches the country code, once both matches fits, it replace the first Year_first previously filled in the species_country file with the year first.
  #Additionally it prints for which species and which countries has the year changed.(Mostly to ensure is working.)
  
  for(i in 1:nrow(exemptions)){
    for (j in 1:nrow(species_country)) {
      try({ 
        if (exemptions$Species_name[i] == species_country$Species_nr[j]){
          if (exemptions$Country_code[i] == species_country$Country_code[j]){

            print(paste0("the species ", exemptions$Sci_name[i], " in ", exemptions$Country_name[i], " has changed from ", species_country$Year_last[j] , " to ", exemptions$Year_last[i]))
            
            species_country$Year_first[j] <- exemptions$Year_first[i]
            species_country$Year_last[j] <- exemptions$Year_last[i]
            
          }
        }
      }, silent = T)
    }
  }
  
  species_country$`Population size (geomean)` = as.numeric(species_country$`Population size (geomean)`)
  species_country$`Population size_min` = as.numeric(species_country$`Population size_min`)
  species_country$`Population size_max` = as.numeric(species_country$`Population size_max`)
  #Now we are checking whether our loops worked properly or not. each loop is writen to give you an error message locating the indices TT folder or 
  #the species and the conflicted years
  #This loop compare the value of Year_first for each species and Year located in first position [1] in its particular indices_TT files
  #The only printed errors should be the pairs "species&countries" within the exemption files.
  
  #for (i in 1:nrow(species_country)) {
  #  try({x<-read.csv(paste0("D:\\WORK\\RSWAN\\RSWAN_tests\\02_2020_PECBMS_data_till2017\\03_EUR\\V5\\working_folder\",species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv"),header = T, sep = ";")
  #  if (species_country$Year_first[i]!= x$Year[1]){
  #print(paste0("Error in", species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv" ))
  #    print(paste0("The correct year for", species_country$SciName[i], " in ", species_country$Country[i], " is ", x$Year[1], " not ", species_country$Year_first[i]))
  
  #     }  
  #  },silent = T)
  #}
  
  #This loop compare the value of Year_first for each species and Year located in first position [nrow(x)] in its particular indices_TT files
  #The only printed errors should be the pairs "species&countries" within the exemption files.
  #for (i in 1:nrow(species_country)) {
  #  try({x<-read.csv(paste0("D:\\WORK\\RSWAN\\RSWAN_tests\\02_2020_PECBMS_data_till2017\\03_EUR\\V5\\working_folder\",species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv"),header = T, sep = ";")
  #    if (species_country$Year_last[i]!= x$Year[nrow(x)]){
  #      #print(paste0("Error in", species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv" ))
  #      print(paste0("The correct year for", species_country$SciName[i], " in ", species_country$Country[i], " is ", x$Year[nrow(x)], " not ", species_country$Year_last[i]))
  #    }
  #  },silent = T)
  #}
  
  
  
  #This export the Species_Country file, directly to the general folder NO THE BACKUP ONE that REMAINS INTACT
  #CHANGE URL!!!!!!!!!!!!!!!!!!!!!!!
  
  write.table(species_country, paste0(general_folder,"Species_Countries.csv"),  
              dec = ".", sep = ";", row.names = FALSE)   
  
  #8310_1_20_indices_TT.csv Netherlands alcedo atis
  #9720_1_4_indices_TT.csv Czech rpublic Galerida cristata
  #10190_1_20_indices_TT.csv Netherlands	Motacilla cinerea
  #11370_1_20_indices_TT.csv Netherlands	Saxicola rubetra
  #11460_1_20_indices_TT.csvNetherlands	Oenanthe oenanthe
  
  # This files, despite been in the combination of the species codes and countries code in the species_working file, do no exist as indices_TT.csv file
  
  # ----------------------------------------------------------------------------
  # END OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  
}

