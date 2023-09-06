#' Download data from the TOV-E database
#'
#' @param dataTable character. The name of the data table to download.
#' @param minYear integer. Earliest year to include in data. 
#' @param maxYear integer. Latest year to include in data. 
#' @param drop_negativeSpp logical. If TRUE (default), species wtih negative and
#' 0 IDs get dropped from the data. If FALSE, all species are kept in data. 
#'
#' @return a tibble containing the requested dataset. 
#' @export
#'
#' @examples

downloadData_TOVE <- function(dataTable, minYear, maxYear, drop_negativeSpp = TRUE){
  
  ## Sort drivers
  sort(unique(odbcListDrivers()[[1]]))
  
  ## Connect to database
  con <- DBI::dbConnect(odbc(),
                        Driver   = "SQL server", 
                        Server   = "ninsql07.nina.no",
                        Database = "TOVTaksering",
                        Trusted_Connection = "True")
  
  ## Download specified data
  TOVE_data <- dplyr::tbl(con, dataTable) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(Site, Species, Year) %>%
    dplyr::filter(Year >= minYear,
                  Year <= maxYear)
  
  ## Optional: remove negative species IDs
  if(drop_negativeSpp){
    TOVE_data <- TOVE_data %>%
      dplyr::filter(Species > 0)
  }
  
  ## Return downloaded data
  return(TOVE_data)
}
