#' List species for direct upload to Nature Index database
#'
#' This function requires the file "NorwegianSppNames.csv" in the "data folder, as well as an open connection to Nature Index database. 
#'
#' @param Spp_selection a tibble containing information on species included in 
#' PECMBS, MSI, and NI. Output of listSpecies_PECBMS().
#' @param Spp_exclude a character vector listing species from the original PECBMS
#' compilation which should not be directly uploaded/updated in NI database.
#'
#' @return a tibble containing information on species that can be updated directly
#' in the Nature Index database based on TRIM analysis results. 
#' @export
#'
#' @examples
#' 
listSpecies_NI <- function(Spp_selection, Spp_exclude){
  
  directNI <- Spp_selection %>%
    dplyr::filter(NI & !(Species %in% Spp_exclude))
  
  
  ## Add Norwegian names for Nature index species
  NorwegianNames <- read.csv("data/NorwegianSppNames.csv")
  directNI <- directNI %>%
    dplyr::left_join(NorwegianNames, by = c("EURINGCode", "Species"))
  

  ## Check for missing species
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
  
  
  ## Retrieve NI indicator IDs for all relevant species
  myIndicators <- NIcalc::getIndicators() %>%
    dplyr::rename(Species_Norwegian_NI = name)
  
  directNI <- directNI %>%
    dplyr::left_join(myIndicators, by = "Species_Norwegian_NI")
  
  
  ## Check that all relevant species have an id in the Nature Index database
  if(any(is.na(directNI$id))){
    NAs <- subset(directNI, is.na(id))
    message(paste0("Indicator ID is missing for the following ", nrow(NAs), " species:"))
    print(NAs$Species)
    
    message("If these are new species to be included in NI, you will have to request their addition to the Nature Index database.")
    message("If these species already exist in the NI database, you may be missing indicator access.")
  }
  
  return(directNI)
}
