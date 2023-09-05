#' List information for species to be included in PECBMS assessment
#'
#' @param file_sppInfo_PECBMS character string. File path for information on
#' PECBMS species. Supported formats: CSV, RDS.
#' @param file_sppList_PECBMS character string. File path for list of included
#' PECBMS species. Supported formats: CSV, RDS.
#'
#' @return a tibble with columns 'EuringCode', 'Year_First', 'TFC_Basetime',
#' 'TCF_Slope_From' and 'Stratum'.
#' @export
#'
#' @examples
#' 
listSpecies_PECBMS <- function(file_sppInfo_PECBMS, file_sppList_PECBMS){
  
  ## Read in species information
  if(grepl(".rds", file_sppInfo_PECBMS)){
    Species_info <- readRDS(file_sppInfo_PECBMS)
  }else{
    
    if(grepl(".csv", file_sppInfo_PECBMS)){
      Species_info <- read.csv2(file_sppInfo_PECBMS)
    }else{
      
      stop("Invalid file format for file_sppInfo_PECBMS. Supported formats: RDS, CSV.")
    }
  }
  
  ## Read in PECBMS species list
  if(grepl(".rds", file_sppInfo_PECBMS)){
    Spp_selection <- readRDS(file_sppList_PECBMS)
  }else{
    
    if(grepl(".csv", file_sppInfo_PECBMS)){
      SSpp_selection <- read.csv2(file_sppList_PECBMS)
    }else{
      
      stop("Invalid file format for file_sppList_PECBMS. Supported formats: RDS, CSV.")
    }
  }
  
  ## Merge and return information
  Spp_selection <- Spp_selection %>%
    tibble::as_tibble() %>%
    dplyr::left_join(Species_info, by = "EURINGCode") %>%
    dplyr::select(EURINGCode, Year_First, TCF_Basetime, TCF_Slope_From) %>%
    dplyr::mutate(Stratum = 63)
  
  return(Spp_selection)
}