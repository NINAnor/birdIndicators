#' List species and double-check information for upload to Nature Index database
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
  
  ## Remove irrelevant species & columns
  sppNI <- Spp_selection %>%
    dplyr::filter(NI & !(Species %in% Spp_exclude)) %>%
    dplyr::select(EURINGCode, Species, indicatorName, indicatorId, dataType,
                  dataUse_direct_SN, dataUse_direct_NN, dataUse_direct_N,
                  dataUse_expert_N) %>%
    dplyr::rename(dataType_NI2020 = dataType)
  
  
  ## Add Norwegian names as per NI database for Nature index species
  NorwegianNames <- read.csv("data/NorwegianSppNames.csv")
  sppNI <- sppNI %>%
    dplyr::left_join(NorwegianNames, by = c("EURINGCode", "Species"))
  

  ## Cross-check Norwegian names in species list with those registered in the database
  
  # In list but missing in NI database
  if(any(is.na(sppNI$Species_Norwegian_NI))){
    NAs <- subset(sppNI, is.na(Species_Norwegian_NI))
    message(paste0("Norwegian species names are missing from the Nature Index database for the following ", nrow(NAs), " species:"))
    print(NAs$Species)
    
    message("If these are new species to be included in NI, please add their Norwegian name to the file data/NorwegianSppNames.csv")
  }
  
  # Different spelling of Norwegian names in list and NI database
  if(any(sppNI$indicatorName != sppNI$Species_Norwegian_NI)){
    
    sub <- subset(sppNI, indicatorName != Species_Norwegian_NI)
    message("The following species have different spelling of Norwegian names in the species list vs. NI database:")
    print(sub[,c("indicatorName", "Species_Norwegian_NI")])
    
    message("This is not a problem. Names as per NI database will be used from here on.")
    
    sppNI <- sppNI %>%
      dplyr::mutate(indicatorName = ifelse(indicatorName != Species_Norwegian_NI, Species_Norwegian_NI, indicatorName))
  }
  
  # Registered in NI database but missing from list
  if(!all(NorwegianNames$Species %in% sppNI$Species)){
    NAs <- subset(NorwegianNames, !(Species %in% sppNI$Species))
    message(paste0("The list of Norwegian names includes", nrow(NAs), " species that are part of Nature Index but for which there seems to be no entry in sppNI:"))
    print(NAs$Species)
    
    message("Please ensure that all necessary species are included in sppNI (check dependencies of function makeSpeciesLists).")
    
  }
  
  
  ## Retrieve NI indicator IDs from database for cross-check
  myIndicators <- NIcalc::getIndicators() %>%
    dplyr::rename(Species_Norwegian_NI = name)
  
  sppNI <- sppNI %>%
    dplyr::left_join(myIndicators, by = "Species_Norwegian_NI")
  
  
  ## Cross-check indicator id's in species list with those registered in the database
  
  # Species with missing indicator id's (not included in subset of indicators user has access to in NI database)
  if(any(is.na(sppNI$id))){
    NAs <- subset(sppNI, is.na(id))
    message(paste0("Indicator ID is missing for the following ", nrow(NAs), " species:"))
    print(NAs$Species)
    
    message("If these are new species to be included in NI, you will have to request their addition to the Nature Index database.")
    message("If these species already exist in the NI database, you may be missing indicator access.")
  }
  
  # Conflicting information on indicator id's in species list and NI database
  if(any(sppNI$indicatorId != sppNI$id)){
    
    sub <- subset(sppNI, indicatorId != id)
    message("The following species have conflicting information on indicator id in the species list vs. NI database:")
    print(sub[,c("indicatorId", "id")])
    
    message("Indicator ids as per NI database will be used from here on.")
    
    sppNI <- sppNI %>%
      dplyr::mutate(indicatorId = ifelse(indicatorId != id, id, indicatorId))
  }
  
  ## Drop duplicate columns
  sppNI <- sppNI %>%
    dplyr::select(-Species_Norwegian_NI, -id)
  
  return(sppNI)
}
