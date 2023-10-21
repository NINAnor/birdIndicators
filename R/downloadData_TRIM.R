#' Download Trim data from the Hekkefuglovervåkingen database
#'
#' @param minYear integer. Earliest year to include in data. 
#' @param maxYear integer. Latest year to include in data. 
#' @param drop_negativeSpp logical. If TRUE (default), species with negative and
#' 0 IDs get dropped from the data. If FALSE, all species are kept in data. 
#' @param DriverName character. Name of the driver to use to access NINA database. 
#' The default ("SQL Server") works for sessions run locally from NINA computers
#' while "FreeTDS" is needed when running from the NINA RStudio server.
#'
#' @return a tibble containing the Trim data for the relevant years including
#' EURING codes.
#' @export
#'
#' @examples

downloadData_TRIM <- function(minYear, maxYear, drop_negativeSpp = TRUE, DriverName = "SQL Server"){
  
  ## Sort drivers
  sort(unique(odbcListDrivers()[[1]]))
  
  ## Connect to database
  con <- DBI::dbConnect(odbc(),
                        Driver   = DriverName, 
                        Server   = "ninsql07.nina.no",
                        Database = "TOVTaksering",
                        port = 1433,
                        Trusted_Connection = "True")
  
  ## Download Trim data
  Trim_data <- dplyr::tbl(con, "TrimResults") %>%
    tibble::as_tibble() %>%
    dplyr::arrange(Site, Species, Year) %>%
    dplyr::filter(Year >= minYear,
                  Year <= maxYear)
  
  ## Download species data
  Spp_data <- dplyr::tbl(con, "Art") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Species = as.numeric(ArtsID)) %>%
    dplyr::select(Species, EURINGCode, Artsnavn_Lat, FK_Kode_Flokk) %>%
    dplyr::arrange(Species) 
  
  ## Optional: remove negative species IDs
  if(drop_negativeSpp){
    Trim_data <- Trim_data %>%
      dplyr::filter(Species > 0)
    Spp_data <- Spp_data %>%
      dplyr::filter(Species > 0)
  }
  
  ## Add the EURING codes to Trim data
  Trim_data_EURING <- Trim_data %>% 
    dplyr::left_join(Spp_data, by = "Species")
  
  ## Count species with missing EURING codes
  noEURING <- Trim_data_EURING %>% 
    dplyr::filter(is.na(EURINGCode)) %>% 
    dplyr::distinct(Artsnavn_Lat) %>%
    nrow()
  
  ## Remove species with missing EURING codes
  Trim_data_EURING <- Trim_data_EURING %>% 
    dplyr::filter(!is.na(EURINGCode)) %>%
    dplyr::mutate(Species_nr = as.numeric(EURINGCode))
  
  ## Remove whitespace from species names
  Trim_data_EURING <- Trim_data_EURING %>%
    dplyr::mutate(Spp_name = stringr::str_trim(Artsnavn_Lat)) %>%
    dplyr::select(-Artsnavn_Lat)
  
  ## Print dataset information
  message(paste0(noEURING, " species were removed as they have no EURING code."))
  message("")
  message(paste0("Number of species in output data: ", length(unique(Trim_data_EURING$Species))))
  message(paste0("Number of sites in output data: ", length(unique(Trim_data_EURING$Site))))
  
  ## Return downloaded data
  return(Trim_data_EURING)
}
