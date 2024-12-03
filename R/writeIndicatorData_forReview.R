#' Write newly calculated NI indicator data to CSV for manual expert review
#'
#' @param dir character. Name of directory into which to save CSV files.
#' @param UpdatedIndicator_data a list containing updated indicator data ready
#' formatted for upload to NI database. 
#'
#' @return one .csv file per species within the specified directory containing
#' updated indicator data ready for review. 
#' @export
#'
#' @examples
writeIndicatorData_forReview <- function(dir, UpdatedIndicator_data){
  
  if(!dir.exists(file.path(dir))){
    dir.create(file.path(dir))
  }
  
  for(i in 1:length(UpdatedIndicator_data)){
    
    if(!is.null(UpdatedIndicator_data[[i]])){
      indID <- UpdatedIndicator_data[[i]]$indicatorValue$indicatorId[1]
      
      readr::write_excel_csv(UpdatedIndicator_data[[i]]$indicatorValues,
                             file = paste0(dir, "/IndID_", indID, ".csv"))
    }
  }
  
}