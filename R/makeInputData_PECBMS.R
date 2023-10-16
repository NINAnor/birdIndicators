#' Set up input data for PECBMS Trim analysis
#'
#' @param Trim_data tibble containing filtered and formatted data from 
#' Hekkefuglovervåkingen. Required columns: Site, Year, Count, EURINGCode, 
#' Species_nr. Output of function downloadData_TRIM(). 
#' @param Spp_selection tibble containing information on species to include in
#' PECBMS analyses. Required columns: EURINGCode, Year_First. Output of function
#' listSpecies_PECBMS(). 
#' @param convertNA logical. If TRUE (default), -1 count values are converted to 
#' NA. If FALSE, they remain as -1. 
#' @param save_allSppData logical. If TRUE (default), saves the data for all
#' species combined as .csv and .rds in working directory.
#' @param returnData logical. If TRUE, returns a nested tibble containing Trim
#' input data for all species that are part of PECBMS assessment. If FALSE 
#' (default), no object is returned. 
#'
#' @return One .csv file per species containing input data for PECBMS assessment
#' as per rtrim-shell format request saved in subfolder "PECBMS_Files". If
#' returnData = TRUE, additionally returns a nested tibble containing the same
#' data.
#' @export
#'
#' @examples

makeInputData_PECBMS <- function(Trim_data, Spp_selection, convertNA = TRUE, save_allSppData = TRUE, returnData = FALSE){
  
  
  ## Selecting the desired species by EURING code
  PECBMS_data <- Trim_data %>%
    dplyr::mutate(EURINGCode = Species_nr) %>%
    dplyr::inner_join(Spp_selection, by = "EURINGCode") 
  
  ## Optional: save complete dataset for PECBMS species
  if(save_allSppData){
    write.csv2(subset(PECBMS_data, PECBMS), paste0("Selected_species_to_PECBMS_", lubridate::year(Sys.Date()), ".csv"), row.names = F)
    saveRDS(subset(PECBMS_data, PECBMS), file = paste0("Selected_species_to_PECBMS_", lubridate::year(Sys.Date()), ".rds"))
  }
  
  ## Optional: convert -1 counts to NA
  if(convertNA){
    PECBMS_data$Count[which(PECBMS_data$Count == -1)] <- NA
  }
  
  ## Filter and reformat input data to match required format
  input_Trim <- PECBMS_data %>% 
    dplyr::select(Site, Year, EURINGCode, Count, Year_First) %>% # Drop unnecessary columns
    dplyr::rename(year = Year,
                  site = Site,
                  count = Count) %>% # Rename columns
    dplyr::mutate(Year_First = ifelse(Year_First == 2006, 2007, Year_First)) %>% # Exclude 2006 as first year
    dplyr::filter(year >= Year_First) %>% # Drop years prior to designated first year
    dplyr::select(site, year, count, EURINGCode) %>% # Re-arrange columns
    dplyr::arrange(EURINGCode, site, year) %>% # Sort by columns
    dplyr::group_by(EURINGCode) %>% # Group by species
    tidyr::nest()
  
  ## Save species-specific input files as per rtrim-shell format request
  for(i in 1:length(input_Trim$EURINGCode)){
    write.csv2(input_Trim$data[[i]], paste0("PECBMS_Files/BMP_", input_Trim$EURINGCode[i],"_1_63_counts.csv"), row.names = FALSE)
  }
  
  ## Optional: return assembled data
  if(returnData){
    return(input_Trim)
  }
  
}