#' List EURING codes and relevant data for species for different targets
#'
#' @param Trim_data a tibble containing the Trim data for the relevant years including
#' EURING codes.
#' @param spp_DataPath character. The path to the folder holding the required 
#' files 'species_table.rds', 'ReportSpeciesPECBMS2020.rds', and
#' 'PECBMS_species_list_2022.rds'.
#' @return a list containing a dataframe of EURING codes with inclusion info
#' for the different workflow targets and relevant information and a list of
#' lists of EURING codes for each target. 
#' @export
#'
#' @examples

makeSpeciesLists <- function(Trim_data, spp_DataPath = "data"){
  
  ## Read in species table
  file_sppTable_ALL <- paste0(spp_DataPath, "/species_table.rds")
  sppTable_ALL <- readRDS(file_sppTable_ALL)
  
  
  ## Load additional information on PECBMS species
  file_sppInfo_PECBMS <- paste0(spp_DataPath, "/ReportSpeciesPECBMS2020.rds")
  file_sppList_PECBMS <- paste0(spp_DataPath, "/PECBMS_species_list_2022.rds")
  
  sppInfo_PECBMS <- listSpecies_PECBMS(file_sppInfo_PECBMS = file_sppInfo_PECBMS, 
                                       file_sppList_PECBMS = file_sppList_PECBMS)

  
  ## Merge PECBMS information into species list
  Spp_selection <- sppTable_ALL %>%
    dplyr::left_join(sppInfo_PECBMS, by = "EURINGCode") %>%
    dplyr::mutate(Stratum = 63)
  
  
  ## Complete Spp_selection for non PECBMS species
  # Extract the first year with data per species from raw data
  sppMinYears <- Trim_data %>%
    dplyr::filter(Count > 0) %>%
    dplyr::group_by(Spp_name) %>%
    dplyr::summarise(minYear = min(Year, na.rm = TRUE)) %>%
    dplyr::rename(Species = Spp_name)
  
  # Add the species-specific starting year ("Year_from" variable)
  Spp_selection <- Spp_selection %>%
    dplyr::left_join(sppMinYears, by = "Species") %>%
    dplyr::mutate(Year_First = ifelse(!is.na(Year_First), Year_First, minYear),
                  TCF_Basetime = 2008, 
                  TCF_Slope_From = ifelse(!is.na(TCF_Slope_From), TCF_Slope_From, minYear)) %>% 
    dplyr::select(-minYear)
  
  # TODO: Check with Diego about criteria for setting Year_First and TCF_Slope_From for non-PECBMS species
  

  ## List EURING codes for species for different outlets
  sppLists <- list(
    ALL = sppTable_ALL$EURINGCode,
    PECBMS = subset(sppTable_ALL, PECBMS)$EURINGCode,
    NI = subset(sppTable_ALL, NI)$EURINGCode,
    MSI = subset(sppTable_ALL, MSI)$EURINGCode,
    MSI_farmland = c(4930, 5410, 9760, 9920, 15820, 18570, 10200),
    MSI_forest = c(10090, 11220, 10840, 12000, 14610, 13140, 14540, 17100, 
                   15390, 8630, 12020, 14860, 8760, 12590, 13120, 12760, 12770,
                   16360, 12010, 13110, 13350, 11870, 14420, 10990),
    MSI_mountain = c(3290, 3300, 4850, 10110, 11060, 11460, 11860, 18470),
    MSI_wetlands = c(5380, 5460, 5480, 5540, 5190, 10170, 18770)
  )
  
  
  ## Assemble and return data and lists
  return(list(sppData = Spp_selection,
              sppLists = sppLists))
}

